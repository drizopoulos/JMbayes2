simulateJoint <- function (alpha = 0.5, Dalpha = 0, n = 500,
                           mean.Cens = 7) {
    # if alpha = 0, mean.Cens = 35
    library("splines")
    library("MASS")
    K <- 15  # number of planned repeated measurements per subject, per outcome
    t.max <- 10 # maximum follow-up time

    ################################################

    # parameters for the linear mixed effects model
    betas <- c("Intercept" = 6.94, "Time1" = 1.30, "Time2" = 1.84, "Time3" = 1.82)
    sigma.y <- 0.6 # measurement error standard deviation

    # parameters for the survival model
    gammas <- c("(Intercept)" = -6.7, "Group" = 0.5)
    #alpha <- 0.5 #0.191
    #Dalpha <- 0.0 # -1.064
    phi <- 2 #2.2
    mean.Cens <- 7

    D <- matrix(0, 4, 4)
    D[lower.tri(D, TRUE)] <- c(0.71, 0.33, 0.07, 1.26, 2.68, 3.81, 4.35, 7.62, 5.4, 8)
    D <- D + t(D)
    diag(D) <- diag(D) * 0.5

    ################################################

    Bkn <- c(0, 9)
    kn <- c(2.1, 3.5)

    # design matrices for the longitudinal measurement model
    times <- c(replicate(n, c(0, 0.5, 1, sort(runif(K - 3, 1, t.max)))))
    group <- rep(0:1, each = n/2)
    DF <- data.frame(time = times)
    X <- model.matrix(~ ns(time, knots = kn, Boundary.knots = Bkn),
                      data = DF)
    Z <- model.matrix(~ ns(time, knots = kn, Boundary.knots = Bkn), data = DF)

    # design matrix for the survival model
    W <- cbind("(Intercept)" = 1, "Group" = group)

    ################################################

    # simulate random effects
    b <- mvrnorm(n, rep(0, nrow(D)), D)

    # simulate longitudinal responses
    id <- rep(1:n, each = K)
    eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ]))
    y <- rnorm(n * K, eta.y, sigma.y)

    # simulate event times
    eta.t <- as.vector(W %*% gammas)
    invS <- function (t, u, i) {
        h <- function (s) {
            NS <- ns(s, knots = kn, Boundary.knots = Bkn)
            DNS <- JMbayes:::dns(s, knots = kn, Boundary.knots = Bkn)
            XX <- cbind(1, NS)
            ZZ <- cbind(1, NS)
            XXd <- DNS
            ZZd <- DNS
            f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
            f2 <- as.vector(XXd %*% betas[2:4] + rowSums(ZZd * b[rep(i, nrow(ZZd)), 2:4]))
            exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha + f2 * Dalpha)
        }
        integrate(h, lower = 0, upper = t)$value + log(u)
    }
    u <- runif(n)
    trueTimes <- numeric(n)
    for (i in 1:n) {
        Up <- 50
        tries <- 5
        Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
        while(inherits(Root, "try-error") && tries > 0) {
            tries <- tries - 1
            Up <- Up + 50
            Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
        }
        trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
    }
    na.ind <- !is.na(trueTimes)
    trueTimes <- trueTimes[na.ind]
    W <- W[na.ind, , drop = FALSE]
    long.na.ind <- rep(na.ind, each = K)
    y <- y[long.na.ind]
    X <- X[long.na.ind, , drop = FALSE]
    Z <- Z[long.na.ind, , drop = FALSE]
    DF <- DF[long.na.ind, , drop = FALSE]
    n <- length(trueTimes)

    Ctimes <- runif(n, 0, 2 * mean.Cens)
    Time <- pmin(trueTimes, Ctimes)
    event <- as.numeric(trueTimes <= Ctimes) # event indicator

    ################################################

    # keep the nonmissing cases, i.e., drop the longitudinal measurements
    # that were taken after the observed event time for each subject.
    ind <- times[long.na.ind] <= rep(Time, each = K)
    y <- y[ind]
    X <- X[ind, , drop = FALSE]
    Z <- Z[ind, , drop = FALSE]
    id <- id[long.na.ind][ind]
    id <- match(id, unique(id))

    dat <- DF[ind, , drop = FALSE]
    dat$id <- id
    dat$y <- y
    dat$Time <- Time[id]
    dat$event <- event[id]
    dat <- dat[c("id", "y", "time", "Time", "event")]
    dat.id <- data.frame(id = unique(dat$id), Time = Time,
                         event = event, group = W[, 2])
    dat$group <- dat.id$group[id]

    summary(tapply(id, id, length))
    n
    mean(event)
    summary(dat.id$Time)
    summary(dat$time)

    # true values for parameters and random effects
    trueValues <- list(betas = betas, tau = 1/sigma.y^2, gammas = gammas,
                       alphas = alpha, Dalphas = Dalpha, sigma.t = phi,
                       inv.D = solve(D), b = b)

    # return list
    list(DF = dat, DF.id = dat.id, trueValues = trueValues)
}

fit_hazard <- function (Data) {
    lmeFit <- lme(y ~ ns(time, k = c(2.1, 3.5), B = c(0, 9)), data = Data$DF,
                  random = list(id = pdDiag(form = ~ ns(time, k = c(2.1, 3.5), B = c(0, 9)))),
                  control = lmeControl(opt = "optim", niterEM = 45))
    coxFit <- coxph(Surv(Time, event) ~ group, data = Data$DF.id)

    JM2 <- jm(coxFit, list(lmeFit), time_var = "time")

    ###########################################################

    test <- JM2


    # parameter values
    betas <- test$initial_values$betas
    b <- test$initial_values$b
    gammas <- test$initial_values$gammas
    bs_gammas <- test$initial_values$bs_gammas
    alphas <- test$initial_values$alphas

    # outcome vectors and design matrices
    n <- test$model_data$n
    idT <- test$model_data$idT
    Time_right <- test$model_data$Time_right
    Time_left <- test$model_data$Time_left
    Time_start <- test$model_data$Time_start
    delta <- test$model_data$delta
    which_event <- test$model_data$which_event
    which_right <- test$model_data$which_right
    which_left <- test$model_data$which_left
    which_interval <- test$model_data$which_interval
    W0_H <- test$model_data$W0_H
    W_H <- test$model_data$W_H
    X_H <- test$model_data$X_H
    Z_H <- test$model_data$Z_H
    U_H <- test$model_data$U_H
    W0_h <- test$model_data$W0_h
    W_h <- test$model_data$W_h
    X_h <- test$model_data$X_h
    Z_h <- test$model_data$Z_h
    U_h <- test$model_data$U_h
    W0_H2 <- test$model_data$W0_H2
    W_H2 <- test$model_data$W_H2
    X_H2 <- test$model_data$X_H2
    Z_H2 <- test$model_data$Z_H2
    U_H2 <- test$model_data$U_H2
    log_Pwk <- test$model_data$log_Pwk
    log_Pwk2 <- test$model_data$log_Pwk2

    control <- test$control
    functional_forms_per_outcome <- test$model_info$fun_forms$functional_forms_per_outcome

    system.time({
        for (m in seq_len(M)) {
            if (m == 1) denominator_surv <- logPC_surv(current_bs_gammas, current_gammas,
                                                       current_alphas, tau_bs_gammas)
            # Update bs_gammas
            for (i in seq_along(current_bs_gammas)) {
                proposed_bs_gammas <- current_bs_gammas
                proposed_bs_gammas[i] <- rnorm(1L, current_bs_gammas[i],
                                               scale_bs_gammas[i])
                numerator_surv <- logPC_surv(proposed_bs_gammas, current_gammas,
                                             current_alphas, tau_bs_gammas)
                log_ratio <- numerator_surv - denominator_surv
                if (is.finite(log_ratio) && min(1, exp(log_ratio)) > runif(1)) {
                    current_bs_gammas <- proposed_bs_gammas
                    denominator_surv <- numerator_surv
                    acceptance_bs_gammas[m, i] <- 1
                }
                if (m > 20) {
                    scale_bs_gammas[i] <-
                        robbins_monro_univ(scale = scale_bs_gammas[i],
                                           acceptance_it = acceptance_bs_gammas[m, i],
                                           it = m, target_acceptance = 0.45)
                }
            }
            post_B_tau_bs_gammas <- prior_B_tau_bs_gammas +
                0.5 * c(crossprod(current_bs_gammas, prior_Tau_bs_gammas) %*%
                            current_bs_gammas)
            tau_bs_gammas <- rgamma(1L, post_A_tau_bs_gammas, post_B_tau_bs_gammas)
            ###
            res_bs_gammas[m, ] <- current_bs_gammas
            res_tau_bs_gammas[m] <- tau_bs_gammas
            res_gammas[m, ] <- current_gammas
            ###
        }
    })
    ###########################
    res_bs_gammas[-seq_len(1000L), ]
}


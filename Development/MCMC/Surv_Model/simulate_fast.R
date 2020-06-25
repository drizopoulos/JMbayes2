library("survival")
library("nlme")
library("GLMMadaptive")
library("splines")
data("pbc2", package = "JM")
data("pbc2.id", package = "JM")
source(file.path(getwd(), "R/jm.R"))
source(file.path(getwd(), "R/help_functions.R"))
source(file.path(getwd(), "Development/jm/R_to_Cpp.R"))
source(file.path(getwd(), "Development/jm/PBC_data.R"))
source(file.path(getwd(), "Development/MCMC/Surv_Model/sample_Surv_Funs.R"))

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
    gammas <- c("(Intercept)" = -9.2, "Group" = 0.5, "Age" = 0.05)
    phi <- 2

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
    age <- runif(n, 30, 70)
    DF <- data.frame(time = times)
    X <- model.matrix(~ ns(time, knots = kn, Boundary.knots = Bkn),
                      data = DF)
    Z <- model.matrix(~ ns(time, knots = kn, Boundary.knots = Bkn), data = DF)

    # design matrix for the survival model
    W <- cbind("(Intercept)" = 1, "Group" = group, "Age" = age)

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
                         event = event, group = W[, 2], age = W[, 3])
    dat$group <- dat.id$group[id]

    #summary(tapply(id, id, length))
    #n
    #mean(event)
    #summary(dat.id$Time)
    #summary(dat$time)

    # true values for parameters and random effects
    trueValues <- list(betas = betas, tau = 1/sigma.y^2, gammas = gammas,
                       alphas = alpha, Dalphas = Dalpha, sigma.t = phi,
                       inv.D = solve(D), b = b)

    # return list
    list(DF = dat, DF.id = dat.id, trueValues = trueValues)
}

fit_hazard <- function (Data, center = FALSE) {
    lmeFit <- lme(y ~ ns(time, k = c(2.1, 3.5), B = c(0, 9)), data = Data$DF,
                  random = list(id = pdDiag(form = ~ ns(time, k = c(2.1, 3.5), B = c(0, 9)))),
                  control = lmeControl(opt = "optim", niterEM = 45))
    coxFit <- coxph(Surv(Time, event) ~ group + age, data = Data$DF.id)

    JM2 <- jm(coxFit, list(lmeFit), time_var = "time")
    test <- JM2

    ###########################################################

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
    # id_H is used to repeat the random effects of each subject GK_k times
    id_H <- lapply(X_H, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
    # this is the linear predictor for the longitudinal outcomes evaluated at the
    # Gauss-Kronrod quadrature points
    eta_H <- linpred_surv(X_H, betas, Z_H, b, id_H)
    # Wlong is the design matrix of all longitudinal outcomes according to the specified
    # functional forms per outcome already multiplied with the interaction terms matrix U
    Wlong_H <- create_Wlong(eta_H, functional_forms_per_outcome, U_H)
    if (length(which_event)) {
        id_h <- lapply(X_h, function (x) seq_len(nrow(x[[1]])))
        eta_h <- linpred_surv(X_h, betas, Z_h, b, id_h)
        Wlong_h <- create_Wlong(eta_h, functional_forms_per_outcome, U_h)
    }
    if (length(which_interval)) {
        id_H2 <- lapply(X_H2, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
        eta_H2 <- linpred_surv(X_H2, betas, Z_H, b, id_H2)
        Wlong_H2 <- create_Wlong(eta_H2, functional_forms_per_outcome, U_H2)
    } else {
        Wlong_H2 <- rep(list(matrix(0.0, length(Time_right), 1)), length(W_H))
    }
    environment(log_density_surv) <- environment()
    environment(logPC_surv) <- environment()
    M <- 5000L
    res_bs_gammas <- acceptance_bs_gammas <- matrix(0.0, M, length(bs_gammas))
    vcov_prop_bs_gammas <- test$vcov_prop$vcov_prop_bs_gammas
    scale_bs_gammas <- rep(0.1, length(bs_gammas))
    prior_mean_bs_gammas <- test$priors$mean_bs_gammas
    prior_Tau_bs_gammas <- test$priors$Tau_bs_gammas
    post_A_tau_bs_gammas <- test$priors$A_tau_bs_gammas +
        0.5 * test$priors$rank_Tau_bs_gammas
    prior_B_tau_bs_gammas <- test$priors$B_tau_bs_gammas
    res_tau_bs_gammas <- numeric(M)
    #
    any_gammas <- test$model_data$any_gammas
    res_gammas <- acceptance_gammas <- matrix(0.0, M, length(gammas))
    vcov_prop_gammas <- test$vcov_prop$vcov_prop_gammas
    scale_gammas <- rep(0.1, length(gammas))
    prior_mean_gammas <- test$priors$mean_gammas
    #
    res_alphas <- acceptance_alphas <- lapply(alphas,
                                              function (a) matrix(0, M, length(a)))
    vcov_prop_alphas <- test$vcov_prop$vcov_prop_alphas
    scale_alphas <- lapply(alphas, function (a) a * 0 + 0.1)
    prior_mean_alphas <- unlist(test$priors$mean_alphas, use.names = FALSE)
    ####
    tau_bs_gammas <- 2
    current_bs_gammas <- jitter(bs_gammas, 80)
    current_gammas <- gammas
    current_alphas <- alphas
    if (center) {
        W_h <- scale(W_h, scale = FALSE)
        W_H <- scale(W_H, scale = FALSE)
        W_H2 <- scale(W_H2, scale = FALSE)
        W_bar <- rbind(attr(W_h, "scaled:center"))
    }
    indFast_H <- id_H[[1]]
    indFast_H <- c(indFast_H[-length(indFast_H)] != indFast_H[-1L], TRUE)
    indFast_H2 <- id_H[[1]]
    indFast_H2 <- c(indFast_H2[-length(indFast_H2)] != indFast_H2[-1L], TRUE)
    W0H_bs_gammas <- W0_H %*% current_bs_gammas
    WH_gammas <- W_H %*% current_gammas
    WlongH_alphas <- Wlong_alphas_fun(Wlong_H, current_alphas)
    if (length(which_event)) {
        W0h_bs_gammas <- W0_h %*% current_bs_gammas
        Wh_gammas <- W_h %*% current_gammas
        Wlongh_alphas <- Wlong_alphas_fun(Wlong_h, current_alphas)
    }
    if (length(which_interval)) {
        W0H2_bs_gammas <- W0_H2 %*% current_bs_gammas
        WH2_gammas <- W_H2 %*% current_gammas
        WlongH2_alphas <- Wlong_alphas_fun(Wlong_H2, current_alphas)
    }
    t0 <- proc.time()
    for (m in seq_len(M)) {
        if (m == 1) {
            denominator_surv <-
                logPC_surv2(current_bs_gammas, current_gammas,
                            current_alphas, tau_bs_gammas,
                            W0H_bs_gammas, WH_gammas, WlongH_alphas,
                            W0h_bs_gammas, Wh_gammas, Wlongh_alphas,
                            W0H2_bs_gammas, WH2_gammas, WlongH2_alphas)
        }
        # Update bs_gammas
        for (i in seq_along(current_bs_gammas)) {
            proposed_bs_gammas <- current_bs_gammas
            proposed_bs_gammas[i] <- rnorm(1L, current_bs_gammas[i],
                                           scale_bs_gammas[i])
            proposed_W0H_bs_gammas <- W0_H %*% proposed_bs_gammas
            if (length(which_event)) {
                proposed_W0h_bs_gammas <- W0_h %*% proposed_bs_gammas
            }
            if (length(which_interval)) {
                proposed_W0H2_bs_gammas <- W0_H2 %*% proposed_bs_gammas
            }
            numerator_surv <-
                logPC_surv2(proposed_bs_gammas, current_gammas,
                            current_alphas, tau_bs_gammas,
                            proposed_W0H_bs_gammas, WH_gammas, WlongH_alphas,
                            proposed_W0h_bs_gammas, Wh_gammas, Wlongh_alphas,
                            proposed_W0H2_bs_gammas, WH2_gammas, WlongH2_alphas)
            log_ratio <- numerator_surv - denominator_surv
            if (is.finite(log_ratio) && min(1, exp(log_ratio)) > runif(1)) {
                current_bs_gammas <- proposed_bs_gammas
                W0H_bs_gammas <- proposed_W0H_bs_gammas
                if (length(which_event)) W0h_bs_gammas <- proposed_W0h_bs_gammas
                if (length(which_interval)) W0H2_bs_gammas <- proposed_W0H2_bs_gammas
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
        if (any_gammas) {
            for (i in seq_along(current_gammas)) {
                proposed_gammas <- current_gammas
                proposed_gammas[i] <- rnorm(1L, current_gammas[i],
                                            scale_gammas[i])
                proposed_WH_gammas <- W_H %*% proposed_gammas
                if (length(which_event)) {
                    proposed_Wh_gammas <- W_h %*% proposed_gammas
                }
                if (length(which_interval)) {
                    proposed_WH2_gammas <- W_H2 %*% proposed_gammas
                }
                numerator_surv <-
                    logPC_surv2(current_bs_gammas, proposed_gammas,
                                current_alphas, tau_bs_gammas,
                                W0H_bs_gammas, proposed_WH_gammas, WlongH_alphas,
                                W0h_bs_gammas, proposed_Wh_gammas, Wlongh_alphas,
                                W0H2_bs_gammas, proposed_WH2_gammas, WlongH2_alphas)
                log_ratio <- numerator_surv - denominator_surv
                if (is.finite(log_ratio) && min(1, exp(log_ratio)) > runif(1)) {
                    current_gammas <- proposed_gammas
                    WH_gammas <- proposed_WH_gammas
                    if (length(which_event)) Wh_gammas <- proposed_Wh_gammas
                    if (length(which_interval)) WH2_gammas <- proposed_WH2_gammas
                    denominator_surv <- numerator_surv
                    acceptance_gammas[m, i] <- 1
                }
                if (m > 20) {
                    scale_gammas[i] <-
                        robbins_monro_univ(scale = scale_gammas[i],
                                           acceptance_it = acceptance_gammas[m, i],
                                           it = m, target_acceptance = 0.45)
                }
            }
        }
        # updates alphas
        for (i in seq_along(current_alphas)) {
            for (j in seq_along(current_alphas[[i]])) {
                proposed_alphas <- current_alphas
                proposed_alphas[[i]][j] <- rnorm(1L, current_alphas[[i]][j],
                                                 scale_alphas[[i]][j])
                proposed_WlongH_alphas <- Wlong_alphas_fun(Wlong_H, proposed_alphas)
                if (length(which_event))
                    proposed_Wlongh_alphas <- Wlong_alphas_fun(Wlong_h, proposed_alphas)
                if (length(which_interval))
                    proposed_WlongH2_alphas <- Wlong_alphas_fun(Wlong_H2, proposed_alphas)
                numerator_surv <-
                    logPC_surv2(current_bs_gammas, current_gammas,
                                proposed_alphas, tau_bs_gammas,
                                W0H_bs_gammas, WH_gammas, proposed_WlongH_alphas,
                                W0h_bs_gammas, Wh_gammas, proposed_Wlongh_alphas,
                                W0H2_bs_gammas, WH2_gammas, proposed_WlongH2_alphas)
                log_ratio <- numerator_surv - denominator_surv
                if (is.finite(log_ratio) && min(1, exp(log_ratio)) > runif(1)) {
                    current_alphas <- proposed_alphas
                    WlongH_alphas <- proposed_WlongH_alphas
                    if (length(which_event)) Wlongh_alphas <- proposed_Wlongh_alphas
                    if (length(which_interval)) WlongH2_alphas <- proposed_WlongH2_alphas
                    denominator_surv <- numerator_surv
                    acceptance_alphas[[i]][m, j] <- 1
                }
                if (m > 20) {
                    scale_alphas[i] <-
                        robbins_monro_univ(scale = scale_alphas[[i]][j],
                                           acceptance_it = acceptance_alphas[[i]][m, j],
                                           it = m, target_acceptance = 0.45)
                }
                res_alphas[[i]][m, j] <- current_alphas[[i]][j]
            }
        }
        ###
        res_bs_gammas[m, ] <- current_bs_gammas
        res_tau_bs_gammas[m] <- tau_bs_gammas
        res_gammas[m, ] <- current_gammas
        ###
    }
    t1 <- proc.time()
    ###########################
    res_bs_gammas <- res_bs_gammas[-seq_len(1000L), ]
    res_gammas <- res_gammas[-seq_len(1000L), , drop = FALSE]
    ttt <- seq(0.0, 12, length.out = 500)
    WW <- splineDesign(test$control$knots, ttt,
                       ord = test$control$Bsplines_degree + 1)
    h0 <- matrix(0.0, nrow(res_bs_gammas), length(ttt))
    for (i in seq_len(nrow(res_bs_gammas))) {
        bs_gammas <- res_bs_gammas[i, ]
        eta <- c(WW %*% bs_gammas)
        if (center && any_gammas) {
            gammas <- res_gammas[i, ]
            eta <- eta - c(W_bar %*% gammas)
        }
        h0[i, ] <- exp(eta)
    }
    list(h0 = colMeans(h0), gammas = colMeans(res_gammas),
         alphas = colMeans(res_alphas[[1]]),
         run_time = t1 - t0)
}

################################################################################
################################################################################


N <- 60
res_h0 <- matrix(0.0, N, 500)
res_gam <- matrix(0.0, N, 2)
res_alph <- matrix(0.0, N, 1)
times <- matrix(0.0, N, 3)
for (j in seq_len(N)) {
    Data_n <- simulateJoint(alpha = 0, mean.Cens = 35)
    fit <- fit_hazard(Data_n, center = TRUE)
    res_h0[j, ] <- fit$h0
    res_gam[j, ] <- fit$gammas
    res_alph[j, ] <- fit$alphas
    times[j, ] <- fit$run_time[1:3]
    print(j)
}

ttt <- seq(0.0, 12, length.out = 500)
plot(x = ttt, y = cbind(colMeans(res_h0)), type = "l",
        lty = c(1), col = 1, xlab = "Time", ylim = c(0, 0.0045),
        ylab = "Baseline Hazard Function")
lines(ttt, exp(Data_n$trueValues$gammas[1] + log(Data_n$trueValues$sigma.t) +
                   (Data_n$trueValues$sigma.t - 1) * log(ttt)), col = "red")



colMeans(res_gam)
Data_n$trueValues$gammas[-1]

colMeans(res_alph)
Data_n$trueValues$alphas

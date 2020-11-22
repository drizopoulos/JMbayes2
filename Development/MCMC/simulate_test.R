library("JMbayes2")
simulateJoint <- function (alpha = 0.5, Dalpha = 0, n = 500,
                           mean.Cens = 7) {
    # if alpha = 0, mean.Cens = 35
    library("splines")
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
    b <- MASS::mvrnorm(n, rep(0, nrow(D)), D)

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
    trueValues <- list(betas = betas, tau = 1/sigma.y^2, gammas = gammas[-1L],
                       alphas = alpha, Dalphas = Dalpha, sigma.t = phi,
                       D = D[lower.tri(D, TRUE)], b = b)

    # return list
    list(DF = dat, DF.id = dat.id, trueValues = trueValues)
}

################################################################################
################################################################################

M <- 300
Data <- simulateJoint()
res_bs_gammas <- matrix(as.numeric(NA), M, 12)
res_gammas <- matrix(as.numeric(NA), M, length(Data$trueValues$gammas))
res_alphas <- matrix(as.numeric(NA), M, length(Data$trueValues$alphas))
res_D <- matrix(as.numeric(NA), M, length(Data$trueValues$D))
res_betas <- matrix(as.numeric(NA), M, length(Data$trueValues$betas))
times <- matrix(as.numeric(NA), M, 3)

for (m in seq_len(M)) {
    try_run <- try({
        Data <- simulateJoint()
        lmeFit <- lme(y ~ ns(time, 3), data = Data$DF,
                      random = list(id = pdDiag(form = ~ ns(time, 3))),
                      control = lmeControl(opt = "optim", niterEM = 45))
        coxFit <- coxph(Surv(Time, event) ~ group + age, data = Data$DF.id)

        obj <- jm(coxFit, list(lmeFit), time_var = "time")
    }, silent = TRUE)
    if (!inherits(try_run, "try-error")) {
        res_bs_gammas[m, ] <- obj$statistics$Mean$bs_gammas
        res_gammas[m, ] <- obj$statistics$Mean$gammas
        res_alphas[m, ] <- obj$statistics$Mean$alphas
        res_D[m, ] <- obj$statistics$Mean$D
        res_betas[m, ] <- obj$statistics$Mean$betas
        times[m, ] <- obj$running_time[1:3L]
    }
    print(m)
}

##########

colMeans(res_gammas, na.rm = TRUE) - Data$trueValues$gammas
colMeans(res_alphas, na.rm = TRUE) - Data$trueValues$alphas
colMeans(res_D, na.rm = TRUE) - Data$trueValues$D
colMeans(res_betas, na.rm = TRUE) - Data$trueValues$betas






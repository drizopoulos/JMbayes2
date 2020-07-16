library("survival")
library("nlme")
library("GLMMadaptive")
library("coda")
library("splines")
library("rbenchmark")
data("pbc2", package = "JM")
data("pbc2.id", package = "JM")
source("./R/jm.R")
source("./R/help_functions.R")
source("./Development/jm/R_to_Cpp.R")
source("./Development/jm/PBC_data.R")
source("./Development/MCMC/Surv_Model/sample_Surv_Funs.R")
Rcpp::sourceCpp('./Development/MCMC/Surv_Model/C++/mcmc.cpp')
#Rcpp::sourceCpp('./Development/MCMC/Surv_Model/C++/mcmc_block.cpp')

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

create_Wlong_mats <- function (model_data, model_info, initial_values, priors,
                               control) {
    betas <- initial_values$betas
    b <- initial_values$b
    gammas <- initial_values$gammas
    bs_gammas <- initial_values$bs_gammas
    alphas <- initial_values$alphas
    # outcome vectors and design matrices
    Time_right <- model_data$Time_right
    Time_left <- model_data$Time_left
    Time_start <- model_data$Time_start
    delta <- model_data$delta
    which_event <- model_data$which_event
    which_right <- model_data$which_right
    which_left <- model_data$which_left
    which_interval <- model_data$which_interval
    W0_H <- model_data$W0_H
    W_H <- model_data$W_H
    X_H <- model_data$X_H
    Z_H <- model_data$Z_H
    U_H <- model_data$U_H
    W0_h <- model_data$W0_h
    W_h <- model_data$W_h
    X_h <- model_data$X_h
    Z_h <- model_data$Z_h
    U_h <- model_data$U_h
    W0_H2 <- model_data$W0_H2
    W_H2 <- model_data$W_H2
    X_H2 <- model_data$X_H2
    Z_H2 <- model_data$Z_H2
    U_H2 <- model_data$U_H2
    # other information
    n <- model_data$n
    idT <- model_data$idT
    log_Pwk <- model_data$log_Pwk
    log_Pwk2 <- model_data$log_Pwk2
    functional_forms_per_outcome <-
        model_info$fun_forms$functional_forms_per_outcome
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
    } else {
        Wlong_h <- rep(list(matrix(0.0, length(Time_right), 1)), length(Wlong_H))
    }
    if (length(which_interval)) {
        id_H2 <- lapply(X_H2, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
        eta_H2 <- linpred_surv(X_H2, betas, Z_H, b, id_H2)
        Wlong_H2 <- create_Wlong(eta_H2, functional_forms_per_outcome, U_H2)
    } else {
        Wlong_H2 <- rep(list(matrix(0.0, length(Time_right), 1)), length(Wlong_H))
    }
    list(Wlong_H = Wlong_H, Wlong_h = Wlong_h, Wlong_H2 = Wlong_H2)
}

Data <- simulateJoint()
lmeFit <- lme(y ~ ns(time, k = c(2.1, 3.5), B = c(0, 9)), data = Data$DF,
              random = list(id = pdDiag(form = ~ ns(time, k = c(2.1, 3.5), B = c(0, 9)))),
              control = lmeControl(opt = "optim", niterEM = 45))
coxFit <- coxph(Surv(Time, event) ~ group + age, data = Data$DF.id)

obj <- test <- jm(coxFit, list(lmeFit), time_var = "time")

model_data <- obj$model_data
model_info <- obj$model_info
initial_values <- obj$initial_values
priors <- obj$priors
vcov_prop <- obj$vcov_prop
control <- obj$control

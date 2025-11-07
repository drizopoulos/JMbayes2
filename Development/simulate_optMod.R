library("JMbayes2")
source("./R/optimal_model.R")

sim_fun1 <- function (n = 100) {
    K <- 15 # number of measurements per subject
    t_max <- 8 # maximum follow-up time
    # we construct a data frame with the design:
    # everyone has a baseline measurement, and then measurements at random
    # follow-up times up to t_max
    DF <- data.frame(id = rep(seq_len(n), each = K),
                     time = rep(seq(0, t_max, length.out = K), n),
                     sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

    # design matrices for the fixed and random effects
    X <- model.matrix(~ sex * time, data = DF)
    Z <- model.matrix(~ time, data = DF)
    betas <- c(-2.2, -0.25, 0.24, -0.05) # fixed effects coefficients
    sigma <- 0.5 # errors' standard deviation
    D11 <- 1.0 # variance of random intercepts
    D22 <- 0.5 # variance of random slopes
    # we simulate random effects
    b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
    # linear predictor
    eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
    # we simulate normal longitudinal data
    DF$y <- rnorm(n * K, mean = eta_y, sd = sigma)

    upp_Cens <- 8 # fixed Type I censoring time
    shape_wb <- 7 # shape Weibull
    alpha <- 0.8 # association coefficients
    gammas <- c("(Intercept)" = -9, "sex" = 0.5)
    W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
    # linear predictor for the survival model
    eta_t <- as.vector(W %*% gammas)
    # to simulate event times we use inverse transform sampling
    # (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want
    # to find t, such that S(t) = u, where S(.) is the survival function, and u a
    # number from the Unif(0, 1) distribution. The function below calculates
    # log(u) - log(S(t)), and for a given u, we want to find t for which it equals
    # zero. We do that below using the uniroot() function
    invS <- function (t, i) {
        # i denotes the subject
        sex_i <- W[i, 2L]
        # h() is the hazard function and we assume a Weibull baseline hazard
        h <- function (s) {
            X_at_s <- cbind(1, sex_i, s, sex_i * s)
            Z_at_s <- cbind(1, s)
            # the linear predictor from the mixed model evaluated at time s
            f <- as.vector(X_at_s %*% betas +
                               rowSums(Z_at_s * b[rep(i, nrow(Z_at_s)), ]))
            exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f * alpha)
        }
        # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
        integrate(h, lower = 0, upper = t)$value + log(u[i])
    }
    # we simulate the event times
    u <- runif(n)
    trueTimes <- numeric(n)
    for (i in seq_len(n)) {
        Up <- 100
        Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
        trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
    }
    # we use fixed Type I right censoring denoting the end of the trial.
    Ctimes <- upp_Cens
    Time <- pmin(trueTimes, Ctimes)
    event <- as.numeric(trueTimes <= Ctimes) # event indicator
    # we keep the longitudinal measurements before the event times
    DF$Time <- Time[DF$id]
    DF$event <- event[DF$id]
    DF$id <- with(DF, ave(id, id, FUN = function (x) sample(1e06, 1) + x))
    DF <- DF[DF$time <= DF$Time, ]
}

sim_fun2 <- function (n = 100) {
    K <- 15 # number of measurements per subject
    t_max <- 8 # maximum follow-up time
    # we construct a data frame with the design:
    # everyone has a baseline measurement, and then measurements at random
    # follow-up times up to t_max
    DF <- data.frame(id = rep(seq_len(n), each = K),
                     time = rep(seq(0, t_max, length.out = K), n),
                     sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

    # design matrices for the fixed and random effects
    X <- model.matrix(~ sex * (time + I(time^2)), data = DF)
    Z <- model.matrix(~ time + I(time^2), data = DF)
    betas <- c(5.2, -0.25, 0.2, 0.1, -0.2, -0.1) # fixed effects coefficients
    sigma <- 0.5 # errors' standard deviation
    D11 <- 0.6 # variance of random intercepts
    D22 <- 0.5 # variance of random slopes
    D33 <- 0.05 # variance of quadratic random slopes
    # we simulate random effects
    b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)),
               rnorm(n, sd = sqrt(D33)))
    # linear predictor
    eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
    # we simulate normal longitudinal data
    DF$y <- rnorm(n * K, mean = eta_y, sd = sigma)

    upp_Cens <- 5.1 # fixed Type I censoring time
    shape_wb <- 7 # shape Weibull
    alpha <- 0.15 # association coefficients
    gammas <- c("(Intercept)" = -12, "sex" = 0.5)
    W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
    # linear predictor for the survival model
    eta_t <- as.vector(W %*% gammas)
    # to simulate event times we use inverse transform sampling
    # (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want
    # to find t, such that S(t) = u, where S(.) is the survival function, and u a
    # number from the Unif(0, 1) distribution. The function below calculates
    # log(u) - log(S(t)), and for a given u, we want to find t for which it equals
    # zero. We do that below using the uniroot() function
    invS <- function (t, i) {
        # i denotes the subject
        sex_i <- W[i, 2L]
        # h() is the hazard function and we assume a Weibull baseline hazard
        h <- function (s) {
            X_at_s <- cbind(1, sex_i, s, s * s, sex_i * s, sex_i * s * s)
            Z_at_s <- cbind(1, s, s * s)
            # the linear predictor from the mixed model evaluated at time s
            f <- as.vector(X_at_s %*% betas +
                               rowSums(Z_at_s * b[rep(i, nrow(Z_at_s)), ]))
            exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f * alpha)
        }
        # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
        integrate(h, lower = 0, upper = t)$value + log(u[i])
    }
    # we simulate the event times
    u <- runif(n)
    trueTimes <- numeric(n)
    for (i in seq_len(n)) {
        Up <- 100
        Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
        trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
    }
    # we use fixed Type I right censoring denoting the end of the trial.
    Ctimes <- upp_Cens
    Time <- pmin(trueTimes, Ctimes)
    event <- as.numeric(trueTimes <= Ctimes) # event indicator
    # we keep the longitudinal measurements before the event times
    DF$Time <- Time[DF$id]
    DF$event <- event[DF$id]
    DF$id <- with(DF, ave(id, id, FUN = function (x) sample(1e06, 1) + x))
    DF <- DF[DF$time <= DF$Time, ]
}

sim_fun3 <- function (n = 100) {
    K <- 15 # number of measurements per subject
    t_max <- 8 # maximum follow-up time
    # we construct a data frame with the design:
    # everyone has a baseline measurement, and then measurements at random
    # follow-up times up to t_max
    DF <- data.frame(id = rep(seq_len(n), each = K),
                     time = rep(seq(0, t_max, length.out = K), n),
                     sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

    # design matrices for the fixed and random effects
    X <- model.matrix(~ time + I(time^2), data = DF)
    Z <- model.matrix(~ time + I(time^2), data = DF)
    betas <- c(5.2, 0.2, 0.1) # fixed effects coefficients
    sigma <- 0.5 # errors' standard deviation
    D11 <- 0.6 # variance of random intercepts
    D22 <- 0.5 # variance of random slopes
    D33 <- 0.05 # variance of quadratic random slopes
    # we simulate random effects
    b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)),
               rnorm(n, sd = sqrt(D33)))
    # linear predictor
    eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
    # we simulate normal longitudinal data
    DF$y <- rnorm(n * K, mean = eta_y, sd = sigma)

    upp_Cens <- 5.1 # fixed Type I censoring time
    shape_wb <- 7 # shape Weibull
    alpha <- 0.15 # association coefficients
    gammas <- c("(Intercept)" = -12, "sex" = 0.5)
    W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
    # linear predictor for the survival model
    eta_t <- as.vector(W %*% gammas)
    # to simulate event times we use inverse transform sampling
    # (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want
    # to find t, such that S(t) = u, where S(.) is the survival function, and u a
    # number from the Unif(0, 1) distribution. The function below calculates
    # log(u) - log(S(t)), and for a given u, we want to find t for which it equals
    # zero. We do that below using the uniroot() function
    invS <- function (t, i) {
        # i denotes the subject
        sex_i <- W[i, 2L]
        # h() is the hazard function and we assume a Weibull baseline hazard
        h <- function (s) {
            X_at_s <- cbind(1, s, s * s)
            Z_at_s <- cbind(1, s, s * s)
            # the linear predictor from the mixed model evaluated at time s
            f <- as.vector(X_at_s %*% betas +
                               rowSums(Z_at_s * b[rep(i, nrow(Z_at_s)), ]))
            exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f * alpha)
        }
        # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
        integrate(h, lower = 0, upper = t)$value + log(u[i])
    }
    # we simulate the event times
    u <- runif(n)
    trueTimes <- numeric(n)
    for (i in seq_len(n)) {
        Up <- 100
        Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
        trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
    }
    # we use fixed Type I right censoring denoting the end of the trial.
    Ctimes <- upp_Cens
    Time <- pmin(trueTimes, Ctimes)
    event <- as.numeric(trueTimes <= Ctimes) # event indicator
    # we keep the longitudinal measurements before the event times
    DF$Time <- Time[DF$id]
    DF$event <- event[DF$id]
    DF$id <- with(DF, ave(id, id, FUN = function (x) sample(1e06, 1) + x))
    DF <- DF[DF$time <= DF$Time, ]
}


################################################################################

training <- do.call('rbind', list(sim_fun1(), sim_fun2(), sim_fun3()))
training$id <- match(training$id, unique(training$id))
testing <- do.call('rbind', list(sim_fun1(), sim_fun2(), sim_fun3()))
testing$id <- match(testing$id, unique(testing$id))
testing2 <- do.call('rbind', list(sim_fun1(), sim_fun2(), sim_fun3()))
testing2$id <- match(testing2$id, unique(testing2$id))


training_id <- training[!duplicated(training$id), ]
CoxFit <- coxph(Surv(Time, event) ~ sex, data = training_id)

fm1 <- lme(y ~ time * sex, data = training, random = ~ time | id)
fm2 <- lme(y ~ sex * ns(time, 3), data = training,
           random = list(id = pdDiag(~ ns(time, 3))))
fm3 <- lme(y ~ poly(time, 2), data = training,
           random = list(id = pdDiag(~ poly(time, 2))))

jointFit1 <- jm(CoxFit, fm1, time_var = "time")
jointFit2 <- jm(CoxFit, fm2, time_var = "time")
jointFit3 <- jm(CoxFit, fm3, time_var = "time")
Models <- list(jointFit1, jointFit2, jointFit3)

T0 <- 3
Data <- testing[ave(testing$time, testing$id, FUN = max) > T0, ]
Data$Time <- T0; Data$event <- 0
Data_after <- Data[Data$time > T0, ]
Preds <- local({
    Data <- Data[Data$time <= T0, ]
    lapply(Models, predict, newdata = Data, newdata2 = Data_after)
})
Preds <- do.call('cbind', lapply(Preds, function (x) x$newdata2$preds[[1L]]))
Obs <- Data_after$y
best_model <- which.min(colMeans((Preds - Obs)^2))

OptModel <- opt_model(Models, Data, T0, cores = 3L)
sapply(OptModel, function (x) sqrt(x$MISE_mod_ave))
mises <- do.call('cbind', lapply(OptModel, function (x) {
    x[["MISEs_ave"]] #+ x[["std_MISEs_vario"]]
}))
best_model_id <- apply(mises, 1L, which.min)
Preds <- do.call('cbind', lapply(OptModel, function (x) x$Preds$newdata2$preds[[1L]]))
Obs <- Data_after$y
id <- match(Data_after$id, unique(Data_after$id))

colMeans((Preds - Obs)^2)[best_model]
mean((Preds[cbind(seq_along(id), best_model_id[id])] - Obs)^2)
weights <- t(apply(1/mises, 1L, function (x) exp(x) / sum(exp(x))))
mean((rowSums(weights[id, ] * Preds) - Obs)^2, na.rm = TRUE)


prs1 <- predict(jointFit1, newdata = Data, return_params_mcmc = TRUE)
prs2 <- predict(jointFit2, newdata = Data, return_params_mcmc = TRUE)
prs3 <- predict(jointFit3, newdata = Data, return_params_mcmc = TRUE)
ppcheck(jointFit1, nsim = 200L, newdata = Data, random_effects = "mcmc",
        params_mcmc = prs1$params_mcmc, type = "average",
        main = c("jointFit1"), pos_legend = c(NA, "left"))
ppcheck(jointFit2, nsim = 200L, newdata = Data, random_effects = "mcmc",
        params_mcmc = prs2$params_mcmc, type = "average",
        main = c("jointFit2"), pos_legend = c(NA, "left"))
ppcheck(jointFit3, nsim = 200L, newdata = Data, random_effects = "mcmc",
        params_mcmc = prs3$params_mcmc, type = "average",
        main = c("jointFit3"), pos_legend = c(NA, "left"))



#################################
#################################
Data2 <- testing2[ave(testing2$time, testing2$id, FUN = max) > T0, ]
Data2_after <- Data2[Data2$time > T0, ]
OptModel <- opt_model(Models, Data2, T0, cores = 3L)
mises <- do.call('cbind', lapply(OptModel, function (x) {
    x[["std_MISEs_ave"]] #+ x[["std_MISEs_vario"]]
}))
best_model_id <- apply(mises, 1L, which.min)
Preds2 <- do.call('cbind', lapply(OptModel, function (x) x$Preds$newdata2$preds[[1L]]))
Obs2 <- Data2_after$y
id2 <- match(Data2_after$id, unique(Data2_after$id))

colMeans((Preds2 - Obs2)^2)[best_model]
mean((Preds2[cbind(seq_along(id2), best_model_id[id2])] - Obs2)^2)
weights <- t(apply(1/mises, 1L, function (x) exp(x) / sum(exp(x))))
mean((rowSums(weights[id2, ] * Preds2) - Obs2)^2, na.rm = TRUE)

round(weights, 3)



xx <- do.call('cbind', lapply(OptModel, function (x) {
    x[["MISEs_ave"]] #+ x[["std_MISEs_vario"]]
}))
round(xx, 3)

sqrt(round(xx, 3))["296", , drop = FALSE]
samplePatient <- Data2[Data2$id == 1, ]
samplePatient <- samplePatient[samplePatient$time <= T0, ]
samplePatient$Time <- T0; samplePatient$event <- 0
prs1 <- predict(jointFit1, newdata = samplePatient, return_params_mcmc = TRUE)
prs2 <- predict(jointFit2, newdata = samplePatient, return_params_mcmc = TRUE)
prs3 <- predict(jointFit3, newdata = samplePatient, return_params_mcmc = TRUE)
ppcheck(jointFit1, nsim = 200L, newdata = samplePatient, random_effects = "mcmc",
        params_mcmc = prs1$params_mcmc, type = "average",
        main = c("jointFit1"), pos_legend = c(NA, "topleft"))
points(samplePatient$time, samplePatient$y, cex = 2)
ppcheck(jointFit2, nsim = 200L, newdata = samplePatient, random_effects = "mcmc",
        params_mcmc = prs2$params_mcmc, type = "average",
        main = c("jointFit2"), pos_legend = c(NA, "topleft"))
ppcheck(jointFit3, nsim = 200L, newdata = samplePatient, random_effects = "mcmc",
        params_mcmc = prs3$params_mcmc, type = "average",
        main = c("jointFit3"), pos_legend = c(NA, "topleft"))


sqrt(round(xx, 3))["13", ]

round(weights, 3)

####

models = Models
t0 <- 2
newdata = Data
parallel = "snow"
cores = 1L

object = Models[[2]]
obs =  sims$outcome[[1]]; reps = sims[[1]]
id = attr(sims$outcome[[1]], "id"); times = attr(sims$outcome[[1]], "times")




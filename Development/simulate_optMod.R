library("JMbayes2")
library("lattice")
library("matrixStats")
source("./R/optimal_model.R")

sim_fun1 <- function (n, model = c("mixed", "joint"), K = 30) {
    model <- match.arg(model)
    t_max <- 7 # maximum follow-up time
    # we construct a data frame with the design:
    # everyone has a baseline measurement, and then measurements at random
    # follow-up times up to t_max
    DF <- data.frame(id = rep(seq_len(n), each = K),
                     time = rep(seq(0, t_max, length.out = K), n),
                     sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

    # design matrices for the fixed and random effects
    X <- model.matrix(~ sex * time, data = DF)
    Z <- model.matrix(~ time, data = DF)
    betas <- 20 * c(-2.2, -0.25, 0.24, -0.05) # fixed effects coefficients
    sigma <- 1 # errors' standard deviation
    D11 <- 2 # variance of random intercepts
    D22 <- 1.5 # variance of random slopes
    # we simulate random effects
    b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
    # linear predictor
    eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
    # we simulate normal longitudinal data
    DF$y <- rnorm(n * K, mean = eta_y, sd = sigma)

    if (model != "mixed") {
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
        DF <- DF[DF$time <= DF$Time, ]
    }
    DF$id <- with(DF, ave(id, id, FUN = function (x) sample(1e06, 1) + x))
    DF
}

sim_fun2 <- function (n, model = c("mixed", "joint"), K = 30) {
    model <- match.arg(model)
    t_max <- 7 # maximum follow-up time
    # we construct a data frame with the design:
    # everyone has a baseline measurement, and then measurements at random
    # follow-up times up to t_max
    DF <- data.frame(id = rep(seq_len(n), each = K),
                     time = rep(seq(0, t_max, length.out = K), n),
                     sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

    # design matrices for the fixed and random effects
    X <- model.matrix(~ sex * ns(time, k = c(3, 5), B = c(0, 7)), data = DF)
    Z <- model.matrix(~ ns(time, k = c(3, 5), B = c(0, 7)), data = DF)
    betas <- 10 * c(5.2, -0.25, 0.2, 0.3, 0.5, -0.2, -0.3, -0.5) # fixed effects coefficients
    sigma <- 1 # errors' standard deviation
    D11 <- 2 # variance of random intercepts
    D22 <- 2 # variance of random slopes
    D33 <- 2 # variance of quadratic random slopes
    D44 <- 2 # variance of quadratic random slopes
    # we simulate random effects
    b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)),
               rnorm(n, sd = sqrt(D33)), rnorm(n, sd = sqrt(D44)))
    # linear predictor
    eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
    # we simulate normal longitudinal data
    DF$y <- rnorm(n * K, mean = eta_y, sd = sigma)

    if (model != "mixed") {
        upp_Cens <- 5.1 # fixed Type I censoring time
        shape_wb <- 7 # shape Weibull
        alpha <- 0.15 # association coefficients
        gammas <- c("(Intercept)" = -8, "sex" = 0.5)
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
        DF <- DF[DF$time <= DF$Time, ]
    }
    DF$id <- with(DF, ave(id, id, FUN = function (x) sample(1e06, 1) + x))
    DF
}

sim_fun3 <- function (n, model = c("mixed", "joint"), K = 30) {
    model <- match.arg(model)
    t_max <- 7 # maximum follow-up time
    # we construct a data frame with the design:
    # everyone has a baseline measurement, and then measurements at random
    # follow-up times up to t_max
    DF <- data.frame(id = rep(seq_len(n), each = K),
                     time = rep(seq(0, t_max, length.out = K), n),
                     sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

    # design matrices for the fixed and random effects
    X <- model.matrix(~ time + I(time^2), data = DF)
    Z <- model.matrix(~ time + I(time^2), data = DF)
    betas <- -10 * c(5.2, 0.2, 0.1) # fixed effects coefficients
    sigma <- 1 # errors' standard deviation
    D11 <- 2 # variance of random intercepts
    D22 <- 2 # variance of random slopes
    D33 <- 2 # variance of quadratic random slopes
    # we simulate random effects
    b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)),
               rnorm(n, sd = sqrt(D33)))
    # linear predictor
    eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
    # we simulate normal longitudinal data
    DF$y <- rnorm(n * K, mean = eta_y, sd = sigma)

    if (model != "mixed") {
        upp_Cens <- 5.1 # fixed Type I censoring time
        shape_wb <- 7 # shape Weibull
        alpha <- 0.15 # association coefficients
        gammas <- c("(Intercept)" = -10, "sex" = 0.5)
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
        DF <- DF[DF$time <= DF$Time, ]
    }
    DF$id <- with(DF, ave(id, id, FUN = function (x) sample(1e06, 1) + x))
    DF
}

create_data <- function (n1, n2, n3, K = 30) {
    DF <- do.call('rbind', list(sim_fun1(n1, K = K), sim_fun2(n2, K = K),
                                sim_fun3(n3, K = K)))
    DF$id <- match(DF$id, unique(DF$id))
    DF
}
fit_models <- function (training) {
    #training_id <- training[!duplicated(training$id), ]
    #CoxFit <- coxph(Surv(Time, event) ~ sex, data = training_id)

    fm1 <- lme(fixed = y ~ time * sex, data = training, random = ~ time | id,
               control = lmeControl(opt = "optim"))
    fm2 <- lme(fixed = y ~ sex * nsk(time, 3), data = training,
               random = list(id = pdDiag(form = ~ nsk(time, 3))),
               control = lmeControl(opt = "optim"))
    fm3 <- lme(fixed = y ~ poly(time, 2) * sex, data = training,
               random = list(id = pdDiag(form = ~ poly(time, 2))),
               control = lmeControl(opt = "optim"))
    fm4 <- lme(fixed = y ~ poly(time, 3), data = training,
               random = list(id = pdDiag(form = ~ poly(time, 3))),
               control = lmeControl(opt = "optim"))

    #jointFit1 <- jm(CoxFit, fm1, time_var = "time")
    #jointFit2 <- jm(CoxFit, fm2, time_var = "time")
    #jointFit3 <- jm(CoxFit, fm3, time_var = "time")
    #jointFit4 <- jm(CoxFit, fm4, time_var = "time")

    list(fm1, fm2, fm3, fm4)
}
best_model_test <- function (testing, T0, Dt) {
    Data <- testing[ave(testing$time, testing$id, FUN = max) > T0, ]
    Data$Time <- T0; Data$event <- 0
    Data_before <- Data[Data$time <= T0, ]
    Data_after <- Data[Data$time > T0 & Data$time <= T0 + Dt, ]
    preds <- lapply(Models, IndvPred_lme, newdata = Data_before, newdata2 = Data_after)
    Preds_after <- do.call('cbind', lapply(preds, function (x) x$predicted_y))
    Obs_after <- Data_after$y
    loss <- function (log_w) {
        log_w <- c(log_w, 0)
        weights <- exp(log_w - logSumExp(log_w))
        mean((rowSums(rep(weights, each = nrow(Preds_after)) * Preds_after) - Obs_after)^2)
    }
    log_w <- c(optim(rep(0, ncol(Preds_after) - 1L), loss, method = "BFGS")$par, 0)
    weights <- exp(log_w - logSumExp(log_w))
    list(best_model = which.min(colMeans((Preds_after - Obs_after)^2)),
         weights = weights)
}
individualized_selection <- function (testing, T0, Dt, best_model, weights) {
    Data <- testing[ave(testing$time, testing$id, FUN = max) > T0, ]
    Data$Time <- T0; Data$event <- 0
    Data_before <- Data[Data$time <= T0, ]
    Data_after <- Data[Data$time > T0 & Data$time <= T0 + Dt, ]
    preds <- lapply(Models, IndvPred_lme, newdata = Data_before,
                    newdata2 = Data_after)
    ####
    Preds <- do.call('cbind', lapply(preds, function (x) x$fitted_y))
    Obs <- Data_before$y
    id <- match(Data_before$id, unique(Data_before$id))
    ####
    Preds_after <- do.call('cbind', lapply(preds, function (x) x$predicted_y))
    Obs_after <- Data_after$y
    id_after <- match(Data_after$id, unique(Data_after$id))
    # oracle predictions (per id best model for the after predictions)
    oracle <-
        apply(rowsum((Preds_after - Obs_after)^2, id_after, reorder = FALSE),
              1L, which.min)
    mse_oracle <- mean((Preds_after[cbind(seq_along(id_after), oracle[id_after])] - Obs_after)^2)
    # MSE best model from testing
    mse_best_model <- colMeans((Preds_after[, best_model, drop = FALSE] - Obs_after)^2)
    # MSE best model per id
    keep <- Data_before$time > 0#T0 - Dt
    mse_id <- rowsum((Preds[keep, ] - Obs[keep])^2, id[keep], reorder = FALSE)
    model_id <- apply(mse_id, 1L, which.min)
    mse_indv <- mean((Preds_after[cbind(seq_along(id_after), model_id[id_after])] - Obs_after)^2)
    # MSE weighted average MSEs per id
    #sigmas <- matrix(sapply(Models, "[[", "sigma"), nrow(Preds), ncol(Preds),
    #                 byrow = TRUE)
    #log_w <- rowsum(dnorm(Obs[keep], Preds[keep, ], sigmas[keep, ], TRUE), id[keep],
    #                reorder = FALSE)
    #weights <- exp(log_w - rowLogSumExps(log_w))
    #mse_indv_w <- mean((rowSums(weights[id_after, ] * Preds_after) - Obs_after)^2)
    #tt <- Data_before$time
    #obs_vr <- variogram(Obs[keep], tt[keep], id[keep])$svar[, "diffs2"]
    #pred_vr <- sapply(Models, fitted_variogram)
    #id_vr <- rep(unique(id[keep]), sapply(c(table(id[keep])),
    #                                      function (ni) ncol(combn(ni, 2))))
    #log_w <- rowsum(-abs(obs_vr - pred_vr), id_vr, reorder = FALSE)
    #weights <- exp(log_w - rowLogSumExps(log_w))
    #mse_indv_w <- mean((rowSums(weights[id_after, ] * Preds_after) - Obs_after)^2)
    sigmas <- matrix(sapply(Models, "[[", "sigma"), nrow(Preds), ncol(Preds),
                     byrow = TRUE)
    resids <- (Obs - Preds) / sigmas
    tt <- Data_before$time
    id_vr <- rep(unique(id[keep]), sapply(c(table(id[keep])),
                                          function (ni) ncol(combn(ni, 2))))
    Coefs <- matrix(0, length(unique(id)), length(Models)) #numeric(length(Models))
    for (k in seq_along(Models)) {
        obs_vr <- variogram(resids[keep, k], tt[keep], id[keep])$svar
        DF_vr <- data.frame(id = id_vr, sqDiff = obs_vr[, 2L], lags = obs_vr[, 1L])
        #lme_vr <- lme(sqDiff ~ lags, data = DF_vr, random = ~ 1 | id,
        #              control = lmeControl(opt = 'optim'))
        #coef(summary(lm(sqDiff ~ lags, data = DF_vr)))
        #plot(sqDiff ~ lags, data = DF_vr)
        #abline(lm(sqDiff ~ lags, data = DF_vr), col = "red", lwd = 2)
        #Coefs[k] <-
        #    log10(anova(lm(sqDiff ~ lags, data = DF_vr),
        #                lm(sqDiff ~ poly(lags, 2), data = DF_vr))$`Pr(>F)`[2L])
        Coefs[, k] <-
            sapply(split(DF_vr, DF_vr$id), function (v) {
                log10(cor.test(x = v$lag, y = v$sqDiff, method = "spearman",
                              exact = FALSE)$p.value)
                })
    }
    weights <- exp(Coefs - rowLogSumExps(Coefs))
    mse_indv_w <- mean((rowSums(weights[id_after, ] * Preds_after) - Obs_after)^2)
    #weights <- exp(Coefs - logSumExp(Coefs))
    #mse_indv_w <- mean((rowSums(rep(weights, each = nrow(Preds_after)) * Preds_after) - Obs_after)^2)
    c(mse_oracle = mse_oracle, mse_best_model = mse_best_model,
      mse_indv = mse_indv, mse_indv_w = mse_indv_w)
}

################################################################################

# training <- create_data(150, 150, 150, K = 20)
# testing <- create_data(450, 2, 2, K = 20)
# testing2 <- create_data(2, 2, 450, K = 20)

######
Times <- seq(1.5, 5.5, 0.5)
Dts <- 1
settings <- expand.grid(T0 = Times, Dt = Dts)
M <- 20
sim_results <- array(NA_real_, c(nrow(settings), 5, M))
dnams <-
    list(paste0("T0=", sprintf("%.1f", settings$T0), ", Dt=",
                sprintf("%.1f", settings$Dt)),
         c("Oracle", "Best_Model", "Best_AIC", "Indv", "Weights_Indv"))
for (m in seq_len(M)) {
    res <- matrix(0, nrow(settings), 5, dimnames = dnams)
    training <- create_data(150, 150, 150, K = 20)
    testing <- create_data(350, 100, 50, K = 20)
    testing2 <- create_data(25, 25, 400, K = 20)
    Models <- fit_models(training)
    aic_best <- which.min(sapply(Models, AIC))
    if (FALSE) {
        i = 1
        T0 = settings$T0[i]
        Dt = 1
        best_model = selected_model[[1]]
        weights = selected_model[[2]]
    }

    for (i in seq_len(nrow(res))) {
        selected_model <- best_model_test(testing, settings$T0[i], settings$Dt[i])
        res[i, ] <- individualized_selection(testing2, settings$T0[i],settings$Dt[i],
                                 c(selected_model[[1]], aic_best),
                                 weights = selected_model[[2]])
    }
    sim_results[, , m] <- res
}

plot_data <- vector("list", M)
model_nams <- c("Oracle", "Best_Model", "Best_AIC", "Indv", "Weights_Indv")
for (m in seq_len(M)) {
    dd <- sim_results[, , m]
    plot_data[[m]] <- data.frame(
        MSE = c(dd),
        T0 = factor(rep(Times, 5), labels = paste("T0 =", Times)),
        Model = factor(rep(model_nams, each = length(Times)),
                       levels = model_nams)
    )
}
plot_data <- do.call('rbind', plot_data)
bwplot(MSE ~ Model | T0, data = plot_data, as.table = TRUE,
       subset = Model %in% c("Best_Model", "Weights_Indv"),
       scales = list(y = list(relation = "free")),
       coef = 1.5, pch = "|", do.out = FALSE, fill = "lightgrey",
       par.strip.text = list(cex = 0.8), ylim = c(0, 25),
       par.settings = list(box.umbrella = list(lty = 1)))


Res <- apply(sim_results, 1:2, median, na.rm = TRUE)
dimnames(Res) <- dimnames(res)
Res

Res <- apply(sim_results, 1:2, sd, na.rm = TRUE)
dimnames(Res) <- dimnames(res)
Res


mse_id_before <-
    rowsum((Preds[keep, ] - Obs[keep])^2, id[keep], reorder = FALSE) /
    c(table(id[keep]))
mse_id_after <-
    rowsum((Preds_after - Obs_after)^2, id_after, reorder = FALSE) /
    c(table(id_after))

cor(cbind(mse_id_before, mse_id_after))[1:4, 5:8]



################################################################################
################################################################################
################################################################################



T0 <- 2
Dt <- 1
Data <- testing[ave(testing$time, testing$id, FUN = max) > T0, ]
Data$Time <- T0; Data$event <- 0
Data_before <- Data[Data$time <= T0, ]
Data_after <- Data[Data$time > T0 & Data$time <= T0 + Dt, ]
preds <- lapply(Models, predict, newdata = Data_before, newdata2 = Data_after)
Preds <- do.call('cbind', lapply(preds, function (x) x$newdata$preds[[1L]]))
Preds <- Preds[Data_before$time > T0 - Dt, ]
Obs <- Data$y[Data$time <= T0 & Data$time > T0 - Dt]
best_model_before <- which.min(colMeans((Preds - Obs)^2))
Preds_after <- do.call('cbind', lapply(preds, function (x) x$newdata2$preds[[1L]]))
Obs_after <- Data_after$y
Data_before <- Data_before[Data_before$time > T0 - Dt, ]
best_model_after <- which.min(colMeans((Preds_after - Obs_after)^2))

oracle <-
    apply(rowsum((Preds_after - Obs_after)^2, Data_after$id, reorder = FALSE),
          1L, which.min)
id <- match(Data_after$id, unique(Data_after$id))
# oracle predictions (per id best model for the after predictions)
mean((Preds_after[cbind(seq_along(id), oracle[id])] - Obs_after)^2)
mse_id <- rowsum((Preds - Obs)^2, Data_before$id, reorder = FALSE) /
    c(table(match(Data_before$id, unique(Data_before$id))))
model_id <- apply(mse_id, 1L, which.min)
id <- match(Data_after$id, unique(Data_after$id))
# MSE after for best model from before measurements
colMeans((Preds_after - Obs_after)^2)[best_model_before]
# MSE after for best model from after measurements
colMeans((Preds_after - Obs_after)^2)[best_model_after]
# MSE best model per id
mean((Preds_after[cbind(seq_along(id), model_id[id])] - Obs_after)^2)
# MSE weighted average MSEs per id
weights <- t(apply(1 / mse_id, 1L, function (x) exp(x) / sum(exp(x))))
mean((rowSums(weights[id, ] * Preds_after) - Obs_after)^2)


prepanel = function (x, y) {
    stats <- boxplot.stats(y)$stats
    list(ylim = range(stats))
}

ni <- c(table(match(Data_before$id, unique(Data_before$id))))
mse_id_before <- rowsum((Preds - Obs)^2, Data_before$id, reorder = FALSE) / ni
ni <- c(table(match(Data_after$id, unique(Data_after$id))))
mse_id_after <- rowsum((Preds_after - Obs_after)^2, Data_after$id, reorder = FALSE) / ni

cor(cbind(mse_id_before, mse_id_after))[1:3, 4:6]

JMbayes::IndvPred_lme(fm1, newdata = Data_before[1:4, ], timeVar = "time", times = 2:8)$predicted_y
JMbayes::IndvPred_lme(fm2, newdata = Data_before[1:4, ], timeVar = "time", times = 2:8)$predicted_y
JMbayes::IndvPred_lme(fm3, newdata = Data_before[1:4, ], timeVar = "time", times = 2:8)$predicted_y


Data_before$Preds1 <- Preds[, 1L]
Data_before$Preds2 <- Preds[, 2L]
Data_before$Preds3 <- Preds[, 3L]
Data_before$Preds4 <- Preds[, 4L]

xyplot(y + Preds1 + Preds2 + Preds3 + Preds4 ~ time | factor(id), data = Data_before,
       type = c("p", "l", "l", "l", "l"), distribute.type = TRUE, cex = 0.8, pch = 8,
       auto.key = TRUE, layout = c(5, 5), subset = id %in% 1:25,
       scales = list(y = list(relation = "free")))



Data_after$Preds1 <- Preds_after[, 1L]
Data_after$Preds2 <- Preds_after[, 2L]
Data_after$Preds3 <- Preds_after[, 3L]
Data_after$Preds4 <- Preds_after[, 4L]


xyplot(y + Preds1 + Preds2 + Preds3 + Preds4 ~ time | factor(id), data = Data_after,
       type = c("p", "l", "l", "l", "l"), distribute.type = TRUE, cex = 0.8, pch = 8,
       auto.key = TRUE, layout = c(5, 5), subset = id %in% 1:25,
       scales = list(y = list(relation = "free")))

Data_before[Data_before$id == 84, ]
Data_after[Data_after$id == 84, ]




fitted_variogram <- function (object) {
    n <- length(unique(id))
    X <- model.matrix(terms(object), Data)
    formYz <- formula(object$modelStruct$reStruct[[1]])
    mfZ <- model.frame(terms(formYz), data = object$data)
    TermsZ <- attr(mfZ, "terms")
    mfZ <- model.frame(TermsZ, data = Data)
    Z <- model.matrix(TermsZ, mfZ)
    sigma <- object$sigma
    D <- lapply(pdMatrix(object$modelStruct$reStruct), "*", sigma^2)[[1]]
    out <- vector("list", n)
    for (j in seq_len(n)) {
        rows <- Data$time > T0 - Dt & Data$time <= T0 & Data$id == j
        X. <- X[rows, ]
        Z. <- Z[rows, ]
        V <- Z. %*% D %*% t(Z.) + diag(sigma, nrow(Z.), nrow(Z.))
        v <- diag(V)
        combs <- combn(ncol(V), 2L)
        vars <- apply(combs, 2L, function (ind) v[ind[1L]] + v[ind[2L]])
        covs <- apply(combs, 2L, function (ind) V[ind[1L], ind[2L]])
        X1 <- X.[combs[1L, ], ]
        X2 <- X.[combs[2L, ], ]
        out[[j]] <- vars - 2 * covs + c((X1 - X2) %*% fixef(object))^2
    }
    0.5 * unlist(out)
}

##
Data <- testing2[ave(testing2$time, testing2$id, FUN = max) > T0, ]
Data$Time <- T0; Data$event <- 0
Data_before <- Data[Data$time <= T0, ]
Data_after <- Data[Data$time > T0, ]
preds <- lapply(Models, predict, newdata = Data_before, newdata2 = Data_after)
Preds <- do.call('cbind', lapply(preds, function (x) x$newdata$preds[[1L]]))
Obs <- Data$y[Data$time <= T0]
ff <- function (x) {
    nx <- length(x)
    if (nx > 1L) rep(c(1, 0), c(1, nx - 1)) else 1
}
ind_first <- as.logical(with(Data_after, ave(time, id, FUN = ff)))
Preds_after <- do.call('cbind', lapply(preds, function (x) x$newdata2$preds[[1L]]))
Obs_after <- Data_after$y

mse_id <- rowsum((Preds - Obs)^2, Data_before$id, reorder = FALSE)
model_id <- apply(mse_id, 1L, which.min)
id <- match(Data_after$id, unique(Data_after$id))
# MSEs after
colMeans((Preds_after - Obs_after)^2)
# MSE after for best model from testing1
colMeans((Preds_after - Obs_after)^2)[best_model_after]
# MSE best model per id
mean((Preds_after[cbind(seq_along(id), model_id[id])] - Obs_after)^2)
# MSE weighted average mses per id
weights <- t(apply(1 / mse_id, 1L, function (x) exp(x) / sum(exp(x))))
mean((rowSums(weights[id, ] * Preds_after) - Obs_after)^2, na.rm = TRUE)



cbind(id, Preds_after, Obs_after, Preds_after[cbind(seq_along(id), model_id[id])])









ff <- function (x) {
    nx <- length(x)
    if (nx > 1L) rep(c(1, 0), c(1, nx - 1)) else 1
}
ind_first <- as.logical(with(Data_after, ave(time, id, FUN = ff)))
Obs <- Data_after$y
colMeans((Preds[ind_first, ] - Obs[ind_first])^2)
best_model <- which.min(colMeans((Preds - Obs)[ind_first, ]^2))

OptModel <- opt_model(Models, Data, T0, cores = 3L)
sapply(OptModel, function (x) sqrt(x$MISE_mod_ave))
mises <- do.call('cbind', lapply(OptModel, function (x) {
    x[["MISEs_ave"]] #+ x[["std_MISEs_vario"]]
}))
best_model_id <- apply(mises, 1L, which.min)
Preds <- do.call('cbind', lapply(OptModel, function (x) x$Preds$newdata2$preds[[1L]]))
Obs <- Data_after$y
id <- match(Data_after$id, unique(Data_after$id))

colMeans((Preds - Obs)[ind_first, ]^2)
colMeans((Preds - Obs)[ind_first, ]^2)[best_model]
mean((Preds[cbind(seq_along(id), best_model_id[id])] - Obs)[ind_first]^2)
weights <- t(apply(1/mises, 1L, function (x) exp(x) / sum(exp(x))))
mean((rowSums(weights[id, ] * Preds) - Obs)[ind_first]^2, na.rm = TRUE)


prs1 <- predict(jointFit1, newdata = Data[Data$time <= T0, ], return_params_mcmc = TRUE)
prs2 <- predict(jointFit2, newdata = Data[Data$time <= T0, ], return_params_mcmc = TRUE)
prs3 <- predict(jointFit3, newdata = Data[Data$time <= T0, ], return_params_mcmc = TRUE)
ppcheck(jointFit1, nsim = 200L, newdata = Data[Data$time <= T0, ], random_effects = "mcmc",
        params_mcmc = prs1$params_mcmc, type = "average",
        main = c("jointFit1"), pos_legend = c(NA, "left"))
ppcheck(jointFit2, nsim = 200L, newdata = Data[Data$time <= T0, ], random_effects = "mcmc",
        params_mcmc = prs2$params_mcmc, type = "average",
        main = c("jointFit2"), pos_legend = c(NA, "left"))
ppcheck(jointFit3, nsim = 200L, newdata = Data[Data$time <= T0, ], random_effects = "mcmc",
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




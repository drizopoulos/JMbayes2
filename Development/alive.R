library("JMbayes2")
library("geepack")

sim_fun <- function (n = 300) {
    K <- 15 # number of measurements per subject
    t_max <- 8 # maximum follow-up time
    # we construct a data frame with the design:
    # everyone has a baseline measurement, and then measurements at random
    # follow-up times up to t_max
    DF <- data.frame(id = rep(seq_len(n), each = K),
                     time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
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
    DF <- DF[DF$time <= DF$Time, ]
    # create dropout
    create_dropout <- function (yy) {
        check <- which(yy < -4)[1L]
        if (!is.na(check) && check < length(yy))
            yy[seq(check + 1, length(yy))] <- NA_real_
        yy
    }
    DF$y2 <- with(DF, ave(y, id, FUN = create_dropout))
    new_event_info <- function (d) {
        d$Time2 <- d$Time
        d$event2 <- d$event
        if (any(na_ind <- is.na(d$y2))) {
            d$Time2 <- max(d$time[!na_ind]) + 0.001
            d$event2 <- 0
        }
        d
    }
    DF_lis <- split(DF, DF$id)
    do.call('rbind', lapply(DF_lis, new_event_info))
}

DF <- sim_fun()
na_ind <- !is.na(DF$y2)
DF_id <- DF[!duplicated(DF$id), ]

# prediction data
pred_dat <- with(DF, data.frame(time = seq(min(time), max(time), len = 20),
                                sex = factor(levels(sex)[1L], levels(sex))))

# sample means
pred_dat$sample_means <- predict(loess(y ~ time, data = DF), newdata = pred_dat)
pred_dat$sample_means2 <- predict(loess(y2 ~ time, data = DF), newdata = pred_dat)


# GEE
gee_fit <- geeglm(y ~ sex * poly(time, 2), data = DF, id = id,
                  family = gaussian(), corstr = "independence")
pred_dat$gee_preds <- predict(gee_fit, newdata = pred_dat)
gee_fit2 <- geeglm(y2 ~ sex * poly(time, 2), data = DF, id = id,
                   family = gaussian(), corstr = "independence")
pred_dat$gee_preds2 <- predict(gee_fit2, newdata = pred_dat)


Cox_fit2 <- coxph(Surv(Time2, event2) ~ sex, data = DF_id)
jm_fit2 <- jm(Cox_fit2, lme_fit2, time_var = "time")
jm_object <- jm_fit2
FF <- function (t, betas, bi, data) {
    sex <- as.numeric(data$sex == "female")
    coefs <- list(alpha = c(1.50226707557132, 2.09986103038462),
                  norm2 = c(1, 1902, 2832.05429710383, 4886.4372564099))
    Poly <- poly(t, degree = 2, coefs = coefs)
    X <- cbind(1, sex, Poly, Poly * sex)
    Z <- cbind(1, Poly)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}

sims <- simulate(jm_object, nsim = 10L, process = "event", Fforms_fun = FF)

lme_object = lme_fit2
dropouts <- DF[DF$event2 == 0 & DF$Time2 < 8, ]
i = 1
dropouts$Time2 <- sims$Times[dropouts$id, i]


# LME
alive_preds_lme <- function (lme_object, correct = FALSE) {
    if (correct) {
        X <- model.matrix(delete.response(terms(lme_object)), DF)
        fits <- fitted(lme_object)
        fits[is.na(fits)] <- predict(lme_object, newdata = DF[is.na(fits), ])
        betas_alive <- solve(crossprod(X), crossprod(X, fits))
    } else {
        X <- model.matrix(terms(lme_object), model.frame(terms(lme_object), DF))
        fits <- fitted(lme_object)
        fits <- fits[!is.na(fits)]
        betas_alive <- solve(crossprod(X), crossprod(X, fits))
    }
    X_pred <- model.matrix(delete.response(terms(lme_object)), pred_dat)
    c(X_pred %*% betas_alive)
}
lme_fit <- lme(y ~ sex * poly(time, 2), data = DF,
               random = list(id = pdDiag(~ poly(time, 2))))
pred_dat$lme_preds <- alive_preds_lme(lme_fit)
lme_fit2 <- lme(y2 ~ sex * poly(time, 2), data = DF, na.action = na.exclude,
                random = list(id = pdDiag(~ poly(time, 2))))
pred_dat$lme_preds2 <- alive_preds_lme(lme_fit2)
pred_dat$lme_preds3 <- alive_preds_lme(lme_fit2, correct = TRUE)


with(pred_dat,
     matplot(time, cbind(sample_means, gee_preds, lme_preds,
                         sample_means2, gee_preds2, lme_preds2, lme_preds3),
             type = "l", lty = c(rep(1:2, each = 3), 1),
             col = c(rep(c("black", "red", "blue"), 2), "magenta"),
             ylab = "Fitted Means", xlab = "Follow-up Time"))
legend("left", c("sample means (no dropout)", "GEE (no dropout)", "LME (no dropout)",
                 "sample means (dropout)", "GEE (dropout)", "LME (dropout)"),
       lty = rep(1:2, each = 3), col = rep(c("black", "red", "blue"), 2), bty = "n")
legend("bottomleft", "corrected LME",lty = 1, col = "magenta", bty = "n")


################################################################################
################################################################################
################################################################################




# JM
Cox_fit <- coxph(Surv(Time, event) ~ sex, data = DF_id)
jm_fit <- jm(Cox_fit, lme_fit, time_var = "time")
DF$fitted_jm <- fitted(jm_fit)[[1]]
X <- jm_fit$model_data$X[[1]]
betas_alive <- solve(crossprod(X), crossprod(X, fitted(jm_fit)[[1]]))
X_pred <- model.matrix(delete.response(terms(jm_fit)[[1]]), pred_dat)
pred_dat$jm_preds <- c(X_pred %*% betas_alive)
##
Cox_fit2 <- coxph(Surv(Time2, event2) ~ sex, data = DF_id)
jm_fit2 <- jm(Cox_fit2, lme_fit2, time_var = "time")
DF$fitted_jm2[!is.na(DF$y2)] <- fitted(jm_fit2)[[1]]
X <- jm_fit2$model_data$X[[1]]
betas_alive2 <- solve(crossprod(X), crossprod(X, fitted(jm_fit2)[[1]]))
X_pred <- model.matrix(delete.response(terms(jm_fit2)[[1]]), pred_dat)
pred_dat$jm_preds2 <- c(X_pred %*% betas_alive2)

with(pred_dat,
     matplot(time, cbind(sample_means, gee_preds, jm_preds,
                         sample_means2, gee_preds2, jm_preds2),
             type = "l", lty = rep(1:2, each = 3),
             col = rep(c("black", "red", "blue"), 2),
             ylab = "Fitted Means", xlab = "Follow-up Time"))




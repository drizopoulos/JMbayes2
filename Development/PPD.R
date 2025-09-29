remotes::install_github("drizopoulos/jmbayes2")
library("JMbayes2")
library("rbenchmark")

# Cox model for the composite event death or transplantation
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
pbc2$status2 <- as.numeric(pbc2$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)

# a linear mixed model for log serum bilirubin
fm1 <- lme(log(serBilir) ~ ns(year, 3) * sex, data = pbc2,
           random = list(id = pdDiag(~ ns(year, 3))))

# the joint model
jointFit1 <- jm(CoxFit, fm1, time_var = "year", save_random_effects = TRUE)


jointFit = jointFit1
# Posterior Predictive Checks - Longitudinal Outcome
ppcheck(jointFit, main = "xcfsf")
ppcheck(jointFit, random_effects = "mcmc")
ppcheck(jointFit, random_effects = "prior", Fforms_fun = FF)


ppcheck(jointFit, type = "vario", CI_loess = TRUE)
ppcheck(jointFit, type = "vario", random_effects = "mcmc")
ppcheck(jointFit, type = "vario", random_effects = "prior", Fforms_fun = FF)

ppcheck(jointFit, type = "varia", CI_loess = TRUE)
ppcheck(jointFit, type = "varia", random_effects = "mcmc")
ppcheck(jointFit, type = "varia", random_effects = "prior", Fforms_fun = FF)

ppcheck(jointFit, type = "avera", CI_loess = TRUE)
ppcheck(jointFit, type = "avera", random_effects = "mcmc")
ppcheck(jointFit, type = "avera", random_effects = "prior", Fforms_fun = FF)


FF <- function (t, betas, bi, data) {
    sex <- as.numeric(data$sex == "female")
    NS <- ns(t, k = c(0.9911, 3.9863), B = c(0, 14.10579))
    X <- cbind(1, NS, sex, NS * sex)
    Z <- cbind(1, NS)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}

ppcheck(jointFit, process = "event", Fforms_fun = FF)
ppcheck(jointFit, process = "event", Fforms_fun = FF,
                   random_effects = "prior")


#############################################################################
#############################################################################

# Cox model for the composite event death or transplantation
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)

# a linear mixed model for prothrombin time
fm1 <- lme(prothrombin ~ ns(year, 2) * sex, data = pbc2,
           random = list(id = pdDiag(~ ns(year, 2))), na.action = na.exclude)

# a Student's-t mixed model for prothrombin time
fm2 <- mixed_model(prothrombin ~ ns(year, 2) * sex, data = pbc2,
                   random = ~ ns(year, 2) || id, family = students.t(df = 4),
                   iter_EM = 0)

# the joint model2
jointFit1 <- jm(CoxFit, fm1, time_var = "year")
jointFit2 <- jm(CoxFit, fm2, time_var = "year")

# Posterior Predictive Checks - Longitudinal Outcome
JMbayes2:::ppcheck(jointFit1)
JMbayes2:::ppcheck(jointFit2)


FF <- function (t, betas, bi, data) {
    sex <- as.numeric(data$sex == "female")
    NS <- ns(t, k = c(2.053444), B = c(0, 14.10579))
    X <- cbind(1, NS, sex, NS * sex)
    Z <- cbind(1, NS)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}

JMbayes2:::ppcheck(jointFit1, process = "event", Fforms_fun = FF)
JMbayes2:::ppcheck(jointFit2, process = "event", Fforms_fun = FF)

JMbayes2:::plot_hazard(jointFit1, tmax = 14)
JMbayes2:::plot_hazard(jointFit2, tmax = 14)

#############################################################################
#############################################################################

# Cox model for the composite event death or transplantation
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
pbc2$status2 <- as.numeric(pbc2$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)

# a linear mixed model for log serum bilirubin
fm <- lme(log(serBilir) ~ ns(year, 3) * sex + age, data = pbc2,
           random = list(id = pdDiag(~ ns(year, 3))))

# the joint model
jointFit <- jm(CoxFit, fm, time_var = "year")

CVdats <- create_folds(pbc2, V = 5, id_var = "id")
fit_model <- function (data) {
    library("JMbayes2")
    # data
    data$event <- as.numeric(data$status != "alive")
    data_id <- data[!duplicated(data$id), ]
    CoxFit <- coxph(Surv(years, status2) ~ sex, data = data_id)
    fm <- lme(log(serBilir) ~ ns(year, k = c(0.9911, 3.9863),
                                 B = c(0, 14.10579)) * sex + age, data = data,
              random = list(id = pdDiag(~ ns(year, k = c(0.9911, 3.9863),
                                             B = c(0, 14.10579)))))
    jointFit <- jm(CoxFit, fm, time_var = "year")
}

system.time({
    cl <- parallel::makeCluster(5L)
    Model_folds <- parallel::parLapply(cl, CVdats$training, fit_model)
    parallel::stopCluster(cl)
})

JMbayes2:::ppcheck(Model_folds, newdata = CVdats$testing,
                   random_effects = "prior", nsim = 100L)



FF <- function (t, betas, bi, data) {
    sex <- as.numeric(data$sex == "female")
    age <- data$age
    NS <- ns(t, k = c(0.9911, 3.9863), B = c(0, 14.10579))
    X <- cbind(1, NS, sex, age, NS * sex)
    Z <- cbind(1, NS)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}

ppcheck(Model_folds, newdata = CVdats$testing,
        Fforms_fun = FF, random_effects = "prior")


prs <- mapply(predict, object = Model_folds, newdata = CVdats$testing,
              SIMPLIFY = FALSE)


ppcheck(Model_folds, newdata = CVdats$testing,
        random_effects = lapply(prs, function (pr) attr(pr, "b_mat")))


dynamic_data <- function (data, t0) {
    data <- data[ave(data$year, data$id, FUN = max) > t0, ]
    data[['years']] <- t0
    data[['event']] <- 0
    list(pre = data[data$year <= t0, ], post = data[data$year > t0, ])
}
testing_data <- lapply(CVdats$testing, dynamic_data, t0 = 6)

prs <- mapply(predict, object = Model_folds,
              newdata = lapply(testing_data, "[[", "pre"),
              newdata2 = lapply(testing_data, "[[", "post"),
              MoreArgs = list(return_mcmc = TRUE, return_newdata = TRUE),
              SIMPLIFY = FALSE)

rep_y <- do.call('rbind', lapply(prs, function (pr) attr(pr$newdata2, "mcmc")[[1L]]))
obs_y <- unlist(lapply(prs, function (pr) log(pr$newdata2$serBilir)))
JMbayes2:::ecdf_compare(rep_y, obs_y, percentiles = c(0.2, 0.9))




#############################################################################
#############################################################################

CoxFit <- coxph(Surv(Time, death) ~ 1, data = aids.id)

fm1 <- lme(CD4 ~ ns(obstime, k = c(6, 10)), data = aids,
           random = list(patient = pdDiag(~ ns(obstime, k = c(6, 10)))))
jointFit1 <- jm(CoxFit, fm1, time_var = "obstime")
##
fm2 <- mixed_model(I(CD4^2) ~ ns(obstime, k = 6), data = aids,
           random = ~ ns(obstime, k = 6) || patient, family = poisson())

jointFit2 <- jm(CoxFit, fm2, time_var = "obstime")
##
fm3 <- mixed_model(I(CD4^2) ~ ns(obstime, k = 6), data = aids,
                   random = ~ ns(obstime, k = 6) || patient,
                   family = negative.binomial())

jointFit3 <- jm(CoxFit, fm3, time_var = "obstime")


# Posterior Predictive Checks - Longitudinal Outcome
ppcheck(jointFit1, nsim = 50L)
ppcheck(jointFit1, nsim = 50L, random_effects = "prior", Fforms_fun = FF)
ppcheck(jointFit1, nsim = 50L, type = "ave")
ppcheck(jointFit1, nsim = 50L, type = "varia")
ppcheck(jointFit1, nsim = 50L, type = "vario")

ppcheck(jointFit2, nsim = 50L)
ppcheck(jointFit2, nsim = 50L, random_effects = "prior", Fforms_fun = FF)
ppcheck(jointFit2, nsim = 50L, type = "ave")
ppcheck(jointFit2, nsim = 50L, type = "varia")
ppcheck(jointFit2, nsim = 50L, type = "vario")

ppcheck(jointFit3, nsim = 50L)
ppcheck(jointFit3, nsim = 50L, random_effects = "prior", Fforms_fun = FF)
ppcheck(jointFit3, nsim = 50L, type = "ave")
ppcheck(jointFit3, nsim = 50L, type = "varia")
ppcheck(jointFit3, nsim = 50L, type = "vario", ylim = c(0, 5000))


# Posterior Predictive Checks - Event Outcome
FF <- function (t, betas, bi, data) {
    NS <- ns(t, k = c(6, 10), B = c(0, 18))
    X <- cbind(1, NS)
    Z <- cbind(1, NS)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}

ppcheck(jointFit1, process = "event", Fforms_fun = FF,
        pos_legend = c("bottomright", "topright"))
ppcheck(jointFit2, process = "event", Fforms_fun = FF)
ppcheck(jointFit3, process = "event", Fforms_fun = FF)


CVdats <- create_folds(aids, V = 5, id_var = "patient")
fit_model <- function (data) {
    library("JMbayes2")
    # data
    data_id <- data[!duplicated(data$patient), ]
    CoxFit <- coxph(Surv(Time, death) ~ 1, data = data_id)
    fm <- mixed_model(I(CD4^2) ~ ns(obstime, k = 6, B = c(0, 18)), data = data,
                      random = ~ ns(obstime, k = 6, B = c(0, 18)) || patient,
                      family = poisson())
    jm(CoxFit, fm, time_var = "obstime")
}

system.time({
    cl <- parallel::makeCluster(5L)
    Model_folds <- parallel::parLapply(cl, CVdats$training, fit_model)
    parallel::stopCluster(cl)
})

JMbayes2:::ppcheck(Model_folds, newdata = CVdats$testing,
                   random_effects = "prior", nsim = 50L, Fforms_fun = FF)

FF <- function (t, betas, bi, data) {
    NS <- ns(t, k = 6, B = c(0, 18))
    X <- cbind(1, NS)
    Z <- cbind(1, NS)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}

object = Model_folds; nsim = 30L; seed = NULL
newdata = CVdats$testing
process = "longitudinal"
random_effects = "prior"
Fforms_fun = FF; include_outcome = FALSE
tol = 0.001; iter = 100L
percentiles = c(0.025, 0.975)


JMbayes2:::ppcheck(Model_folds, newdata = CVdats$testing, process = "event",
                   Fforms_fun = FF, random_effects = "prior")



#############################################################################
#############################################################################

CoxFit <- coxph(Surv(Time, death) ~ 1, data = prothros)

fm <- lme(pro ~ ns(time, 3), data = prothro,
           random = list(id = pdDiag(~ ns(time, 3))))

jointFit <- jm(CoxFit, fm, time_var = "time", save_random_effects = TRUE)

# Posterior Predictive Checks - Longitudinal Outcome
JMbayes2:::ppcheck(jointFit)
JMbayes2:::ppcheck(jointFit, random_effects = "mc")
JMbayes2:::ppcheck(jointFit, random_effects = "prior", Fforms_fun = FF)


# Posterior Predictive Checks - Event Outcome
FF <- function (t, betas, bi, data) {
    NS <- ns(t, k = c(0.49283, 2.1547), B = c(0, 11.1078))
    X <- cbind(1, NS)
    Z <- cbind(1, NS)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}

JMbayes2:::ppcheck(jointFit, process = "event", Fforms_fun = FF)


plot(ecdf(prothro$pro), verticals = T)

################################################################################
################################################################################
################################################################################

object = jointFit
nsim = 50L
percentiles = c(0.025, 0.975)

out <- simulate(object, nsim = nsim, process = "longitudinal")
n_outcomes <- length(object$model_data$y)
index <- seq_len(n_outcomes)
for (j in index) {
    y <- object$model_data$y[[j]]
    r1 <- quantile(y, probs = percentiles[1L], na.rm = TRUE)
    r2 <- quantile(y, probs = percentiles[2L], na.rm = TRUE)
    x_vals <- seq(r1, r2, length.out = 500)
    rep_y <- sapply(out, function (x, x_vals) ecdf(x[[j]])(x_vals),
                    x_vals = x_vals)
    F0 <- ecdf(y)(x_vals)
    matplot(x_vals, rep_y, type = "l", lty = 1, col = "lightgrey",
            xlab = object$model_info$var_names$respVars_form[[j]],
            ylab = "Empirical CDF")
    lines(x_vals, F0)
    legend("bottomright", c("replicated data", "observed data"),
           lty = 1, col = c("lightgrey", "black"), bty = "n", cex = 0.8)

    trapezoid_rule <- function (f, x) {
        sum(0.5 * diff(x) * (f[-length(x)] + f[-1]))
    }

    MISE <- mean(apply((rep_y - F0)^2, 2L, trapezoid_rule, x = x_vals))
}










object <- jointFit
nsim = 30
Bspline_dgr <- jointFit$control$Bsplines_degree
knots <- jointFit$control$knots[[1]]
W <- object$model_data$W_h
Times <- object$model_data$Time_right
event <- jointFit$model_data$delta
n <- object$model_data$n

# Parameters
mcmc_betas <- do.call('rbind', object$mcmc$betas1)
mcmc_bs_gammas <- do.call('rbind', object$mcmc$bs_gammas)
has_gammas <- !is.null(object$mcmc$gammas)
if (has_gammas) mcmc_gammas <- do.call('rbind', object$mcmc$gammas)
mcmc_W_std_gammas <- do.call('rbind', object$mcmc$W_std_gammas)
mcmc_alphas <- do.call('rbind', object$mcmc$alphas)
mcmc_Wlong_std_alphas <- do.call('rbind', object$mcmc$Wlong_std_alphas)
mcmc_b <- abind::abind(object$mcmc$b)
b <- ranef(object)
nRE <- ncol(b)
get_D <- function (x, n) {
    m <- matrix(0.0, n, n)
    m[lower.tri(m, TRUE)] <- x
    m <- m + t(m)
    diag(m) <- diag(m) * 0.5
    m
}
xx <- do.call('rbind', object$mcmc$D)
mcmc_D <- array(0.0, c(nRE, nRE, nrow(xx)))
for (m in seq_len(nrow(mcmc_betas))) {
    mcmc_D[, , m] <- get_D(xx[m, ], nRE)
}
opt <- 2

hazard <- function (t, subj, rescale_factor, log_u) {
    W0 <- splineDesign(knots = knots, x = t, ord = Bspline_dgr + 1,
                       outer.ok = TRUE)
    log_h0 <- c(W0 %*% bs_gammas) - rescale_factor
    covariates <- c(W[subj, , drop = FALSE] %*% gammas)
    sex <- as.numeric(pbc2.id$sex[subj] == "female")
    NS <- ns(t, k = c(0.9911, 3.9863), B = c(0, 14.10579))
    X <- cbind(1, NS, sex, NS * sex)
    Z <- cbind(1, NS)
    b_subj <- b[subj, ]
    eta <- c(X %*% betas) + rowSums(Z * rep(b_subj, each = nrow(Z)))
    Wlong <- cbind(eta)
    long <- c(Wlong %*% alphas)
    exp(log_h0 + covariates + long)
}

invS <- function (t, log_u, subj, rescale_factor) {
    integrate(hazard, lower = 0, upper = t, subj = subj,
              rescale_factor = rescale_factor)$value + log_u
}

log_uu <- log(runif(1L))
sbj <- 20
tt <- 1.5
t1 <- tt - 1e-03
t2 <- tt + 1e-03
f1 <- invS(t = t1, log_u = log_uu, subj = sbj, rescale_factor = rescale_factor)
f2 <- invS(t = t2, log_u = log_uu, subj = sbj, rescale_factor = rescale_factor)
(f2 - f1) / (t2 - t1)
hazard(tt, sbj, rescale_factor)

fn <- invS(tt, log_uu, sbj, rescale_factor)
gr <- hazard(tt, sbj, rescale_factor)
tt - fn / gr

nr_root <- function (interval, fn, gr, ..., tol = 0.001, iter = 30L) {
    Low <- low <- interval[1L]
    Upp <- upp <- interval[2L]
    fn_Low <- fn(Low, ...)
    fn_Upp <- fn(Upp, ...)
    if (fn_Upp < 0) return(Upp)
    tt <- tt_old <- sum(interval) / 2
    for (i in seq_len(iter)) {
        ffn <- fn(tt, ...)
        # check convergence
        if (abs(ffn) < tol) return(tt)
        # if proposed value outside interval do bisection, else Newton-Raphson
        ggr <- gr(tt, ...)
        tt <- tt - ffn / ggr
        if (tt < Low || tt > Upp) {
            if (ffn < 0 && ffn > fn_Low) {
                low <- tt_old
                fn_Low <- ffn
            } else {
                upp <- tt_old
            }
            tt <- tt_old <- (low + upp) / 2
        }
    }
    tt
}

nr_root(c(0, Up), invS, hazard, log_u = log_uu, subj = sbj,
        rescale_factor = rescale_factor)
uniroot(invS, interval = c(0, Up), log_u = log_uu, subj = sbj,
        rescale_factor = rescale_factor, tol = 0.001)$root

benchmark(nr = nr_root(c(0, Up), invS, hazard, log_u = log_uu, subj = sbj,
                       rescale_factor = rescale_factor),
          uniroot = uniroot(invS, interval = c(0, Up), log_u = log_uu, subj = sbj,
                            rescale_factor = rescale_factor, tol = 0.001)$root)


valT <- eventT <- matrix(0.0, n, nsim)
indices <- sample(nrow(mcmc_betas), nsim)
for (j in seq_len(nsim)) {
    jj <- indices[j]
    betas <- mcmc_betas[jj, ]
    b <- switch(opt,
                "1" = mcmc_b[, , jj],
                "2" = ranef(object),
                "3" = MASS::mvrnorm(n, mu = rep(0, nRE), mcmc_D[, , jj]))
    bs_gammas <- mcmc_bs_gammas[jj, ]
    if (has_gammas) gammas <- mcmc_gammas[jj, ]
    alphas <- mcmc_alphas[jj, ]
    Wlong_std_alphas <- mcmc_Wlong_std_alphas[jj, ]
    W_std_gammas <- mcmc_W_std_gammas[jj, ]
    rescale_factor <- Wlong_std_alphas + W_std_gammas
    Up <- max(Times) * 1.05
    rep_eventTimes1 <- rep_eventTimes2 <- Times
    for (i in seq_len(n)) {
        log_u <- log(runif(1))
        Root <-
            try(uniroot(invS, interval = c(0, Up), log_u = log_u, subj = i,
                        rescale_factor = rescale_factor)$root,
                silent = TRUE)
        rep_eventTimes1[i] <- if (!inherits(Root, "try-error")) Root else Up
        rep_eventTimes2[i] <- nr_root(c(0, Up), invS, hazard, log_u = log_u,
                                      subj = i,
                                      rescale_factor = rescale_factor)
    }
    valT[, j] <- rep_eventTimes
    eventT[, j] <- as.numeric(rep_eventTimes < Up)
}

cbind(rep_eventTimes1, rep_eventTimes2)


xvals <- seq(0, max(Times), length.out = 500)
plot(range(Times), c(0, 1), type = "n", xlab = "Time",
     ylab = "Cumulative Incidence")
for (k in seq_len(ncol(valT))) {
    #lines(xvals, ecdf(valT[, k])(xvals), col = "lightgrey")
    lines(survfit(Surv(valT[, k], eventT[, k]) ~ 1), col = "lightgrey",
          conf.int = FALSE, fun = "event")
}
lines(survfit(Surv(Times, event) ~ 1), fun = "event")
legend("bottomright", c("replicated data", "observed data"), lty = 1,
       col = c("lightgrey", "black"), bty = "n", cex = 0.8)

JMbayes2:::plot_hazard(jointFit, tmax = 14)

###

object <- jointFit
nsim = 20
Bspline_dgr <- jointFit$control$Bsplines_degree
knots <- jointFit$control$knots[[1]]
W <- object$model_data$W_h
Times <- object$model_data$Time_right
event <- jointFit$model_data$delta
n <- object$model_data$n

# Parameters
mcmc_betas <- do.call('rbind', object$mcmc$betas1)
mcmc_bs_gammas <- do.call('rbind', object$mcmc$bs_gammas)
has_gammas <- !is.null(object$mcmc$gammas)
if (has_gammas) mcmc_gammas <- do.call('rbind', object$mcmc$gammas)
mcmc_W_std_gammas <- do.call('rbind', object$mcmc$W_std_gammas)
mcmc_alphas <- do.call('rbind', object$mcmc$alphas)
mcmc_Wlong_std_alphas <- do.call('rbind', object$mcmc$Wlong_std_alphas)
mcmc_b <- abind::abind(object$mcmc$b)
b <- ranef(object)
nRE <- ncol(b)
get_D <- function (x, n) {
    m <- matrix(0.0, n, n)
    m[lower.tri(m, TRUE)] <- x
    m <- m + t(m)
    diag(m) <- diag(m) * 0.5
    m
}
xx <- do.call('rbind', object$mcmc$D)
mcmc_D <- array(0.0, c(nRE, nRE, nrow(xx)))
for (m in seq_len(nrow(mcmc_betas))) {
    mcmc_D[, , m] <- get_D(xx[m, ], nRE)
}
opt <- 2

invS <- function (t, u, subj, rescale_factor) {
    hazard <- function (t, subj) {
        W0 <- splineDesign(knots = knots, x = t, ord = Bspline_dgr + 1,
                           outer.ok = TRUE)
        log_h0 <- c(W0 %*% bs_gammas) - rescale_factor
        covariates <- if (has_gammas) c(W[subj, , drop = FALSE] %*% gammas) else 0.0
        NS <- ns(t, k = 6, B = c(0, 18))
        X <- cbind(1, NS)
        Z <- cbind(1, NS)
        b_subj <- b[subj, ]
        eta <- c(X %*% betas) + rowSums(Z * rep(b_subj, each = nrow(Z)))
        Wlong <- cbind(eta)
        long <- c(Wlong %*% alphas)
        exp(log_h0 + covariates + long)
    }
    integrate(hazard, lower = 0, upper = t, subj = subj)$value + log(u)
}

valT <- eventT <- matrix(0.0, n, nsim)
indices <- sample(nrow(mcmc_betas), nsim)
for (j in seq_len(nsim)) {
    jj <- indices[j]
    betas <- mcmc_betas[jj, ]
    b <- switch(opt,
                "1" = mcmc_b[, , jj],
                "2" = ranef(object),
                "3" = MASS::mvrnorm(n, mu = rep(0, nRE), mcmc_D[, , jj]))
    bs_gammas <- mcmc_bs_gammas[jj, ]
    if (has_gammas) gammas <- mcmc_gammas[jj, ]
    alphas <- mcmc_alphas[jj, ]
    Wlong_std_alphas <- mcmc_Wlong_std_alphas[jj, ]
    W_std_gammas <- mcmc_W_std_gammas[jj, ]
    rescale_factor <- Wlong_std_alphas + W_std_gammas
    Up <- max(Times) * 1.05
    rep_eventTimes <- Times
    for (i in seq_len(n)) {
        u <- runif(1L, 0.0, 1.0)
        Root <-
            try(uniroot(invS, interval = c(0, Up), u = u, subj = i,
                        rescale_factor = rescale_factor)$root,
                silent = TRUE)
        rep_eventTimes[i] <- if (!inherits(Root, "try-error")) Root else Up
    }
    valT[, j] <- rep_eventTimes
    eventT[, j] <- as.numeric(rep_eventTimes < max(Times)) #* event
}

xvals <- seq(0, max(Times), length.out = 500)
plot(range(Times), c(0, 1), type = "n", xlab = "Time",
     ylab = "Cumulative Incidence")
for (k in seq_len(ncol(valT))) {
    lines(xvals, ecdf(valT[, k])(xvals), col = "lightgrey")
    #lines(survfit(Surv(valT[, k], eventT[, k]) ~ 1), col = "lightgrey",
    #      conf.int = FALSE, fun = "event")
}
lines(survfit(Surv(Times, event) ~ 1), fun = "event")
legend("topleft", c("replicated data", "observed data"), lty = 1,
       col = c("lightgrey", "black"), bty = "n", cex = 0.8)


JMbayes2:::plot_hazard(jointFit)


###

object <- jointFit
nsim = 50
Bspline_dgr <- jointFit$control$Bsplines_degree
knots <- jointFit$control$knots[[1]]
W <- object$model_data$W_h
Times <- object$model_data$Time_right
event <- jointFit$model_data$delta
n <- object$model_data$n

# Parameters
mcmc_betas <- do.call('rbind', object$mcmc$betas1)
mcmc_bs_gammas <- do.call('rbind', object$mcmc$bs_gammas)
has_gammas <- !is.null(object$mcmc$gammas)
if (has_gammas) mcmc_gammas <- do.call('rbind', object$mcmc$gammas)
mcmc_W_std_gammas <- do.call('rbind', object$mcmc$W_std_gammas)
mcmc_alphas <- do.call('rbind', object$mcmc$alphas)
mcmc_Wlong_std_alphas <- do.call('rbind', object$mcmc$Wlong_std_alphas)
mcmc_b <- abind::abind(object$mcmc$b)
b <- ranef(object)
nRE <- ncol(b)
get_D <- function (x, n) {
    m <- matrix(0.0, n, n)
    m[lower.tri(m, TRUE)] <- x
    m <- m + t(m)
    diag(m) <- diag(m) * 0.5
    m
}
xx <- do.call('rbind', object$mcmc$D)
mcmc_D <- array(0.0, c(nRE, nRE, nrow(xx)))
for (m in seq_len(nrow(mcmc_betas))) {
    mcmc_D[, , m] <- get_D(xx[m, ], nRE)
}
opt <- 2

invS <- function (t, u, subj, rescale_factor) {
    hazard <- function (t, subj, rescale_factor) {
        W0 <- splineDesign(knots = knots, x = t, ord = Bspline_dgr + 1,
                           outer.ok = TRUE)
        log_h0 <- c(W0 %*% bs_gammas) - rescale_factor
        covariates <- if (has_gammas) c(W[subj, , drop = FALSE] %*% gammas) else 0.0
        NS <- ns(t, k = c(0.49283, 2.1547), B = c(0, 11.1078))
        X <- cbind(1, NS)
        Z <- cbind(1, NS)
        b_subj <- b[subj, ]
        eta <- c(X %*% betas) + rowSums(Z * rep(b_subj, each = nrow(Z)))
        Wlong <- cbind(eta)
        long <- c(Wlong %*% alphas)
        exp(log_h0 + covariates + long)
    }
    integrate(hazard, lower = 0, upper = t, subj = subj,
              rescale_factor = rescale_factor)$value + log(u)
}

valT <- eventT <- matrix(0.0, n, nsim)
indices <- sample(nrow(mcmc_betas), nsim)
for (j in seq_len(nsim)) {
    jj <- indices[j]
    betas <- mcmc_betas[jj, ]
    b <- switch(opt,
                "1" = mcmc_b[, , jj],
                "2" = ranef(object),
                "3" = MASS::mvrnorm(n, mu = rep(0, nRE), mcmc_D[, , jj]))
    bs_gammas <- mcmc_bs_gammas[jj, ]
    if (has_gammas) gammas <- mcmc_gammas[jj, ]
    alphas <- mcmc_alphas[jj, ]
    Wlong_std_alphas <- mcmc_Wlong_std_alphas[jj, ]
    W_std_gammas <- mcmc_W_std_gammas[jj, ]
    rescale_factor <- Wlong_std_alphas + W_std_gammas
    Up <- max(Times) * 1.05
    rep_eventTimes <- Times
    for (i in seq_len(n)) {
        u <- runif(1L, 0.0, 1.0)
        Root <-
            try(uniroot(invS, interval = c(0, Up), u = u, subj = i,
                        rescale_factor = rescale_factor)$root,
                silent = TRUE)
        rep_eventTimes[i] <- if (!inherits(Root, "try-error")) Root else Up
    }
    valT[, j] <- rep_eventTimes
    eventT[, j] <- as.numeric(rep_eventTimes < max(Times)) #* event
}

xvals <- seq(0, max(Times), length.out = 500)
plot(range(Times), c(0, 1), type = "n", xlab = "Time",
     ylab = "Cumulative Incidence")
for (k in seq_len(ncol(valT))) {
    #lines(xvals, ecdf(valT[, k])(xvals), col = "lightgrey")
    lines(survfit(Surv(valT[, k], eventT[, k]) ~ 1), col = "lightgrey",
          conf.int = FALSE, fun = "event")
}
lines(survfit(Surv(Times, event) ~ 1), fun = "event")
legend("topleft", c("replicated data", "observed data"), lty = 1,
       col = c("lightgrey", "black"), bty = "n", cex = 0.8)


JMbayes2:::plot_hazard(jointFit, tmax = 10)


y <- jointFit$model_data$y[[1]]
x_vals <- seq(min(y), max(y), length.out = 500)
rep_y <- sapply(out, function (x, x_vals) ecdf(x[[1]])(x_vals), x_vals = x_vals)
matplot(x_vals, rep_y, type = "l", lty = 1, col = "lightgrey",
        xlab = "Response Variable", ylab = "Empirical CDF")
lines(x_vals, ecdf(y)(x_vals))
legend("bottomright", c("replicated data", "observed data"), lty = 1,
       col = c("lightgrey", "black"), bty = "n", cex = 0.8)


# Posterior Predictive Checks - Event Outcome
FF <- function (t, betas, bi, data) {
    sex <- as.numeric(data$sex == "female")
    NS <- ns(t, k = c(0.9911, 3.9863), B = c(0, 14.10579))
    X <- cbind(1, NS, sex, NS * sex)
    Z <- cbind(1, NS)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}
FF. <- function (t, betas, bi, data) {
    sex <- as.numeric(data$sex == "female")
    NS <- ns(t, k = c(0.9911, 3.9863), B = c(0, 14.10579))
    X <- cbind(1, NS, sex, NS * sex)
    Z <- cbind(1, NS)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * rep(bi, each = nrow(Z)))
    cbind(eta)
}


Times <- jointFit$model_data$Time_right
event <- jointFit$model_data$delta
out1 <- simulate(jointFit, nsim = 10L, process = "e", Fforms_fun = FF)
out2 <- JMbayes2:::simulate2.jm(jointFit, nsim = 10L, process = "e", Fforms_fun = FF)

system.time(out1 <- simulate(jointFit, nsim = 10L, process = "e", Fforms_fun = FF))
system.time(out2 <- JMbayes2:::simulate1.jm(jointFit, nsim = 10L, process = "e", Fforms_fun = FF.))
system.time(out3 <- JMbayes2:::simulate2.jm(jointFit, nsim = 10L, process = "e", Fforms_fun = FF))

object = jointFit
nsim = 10L
seed = 123L
process = "event"
random_effects = "prior"
Fforms_fun = FF; tol = 0.001; iter = 100L
newdata = NULL
percentiles = c(0.025, 0.975)



xvals <- seq(0, max(Times), length.out = 500)
plot(range(Times), c(0, 1), type = "n", xlab = "Time",
     ylab = "Cumulative Incidence")
for (k in seq_len(ncol(out1$Times))) {
    lines(survfit(Surv(out1$Times[, k], out1$event[, k]) ~ 1), col = "lightgrey",
          conf.int = FALSE, fun = "event")
}
lines(survfit(Surv(Times, event) ~ 1), fun = "event")
legend("bottomright", c("replicated data", "observed data"), lty = 1,
       col = c("lightgrey", "black"), bty = "n", cex = 0.8)

xvals <- seq(0, max(Times), length.out = 500)
plot(range(Times), c(0, 1), type = "n", xlab = "Time",
     ylab = "Cumulative Incidence")
for (k in seq_len(ncol(out2$Times))) {
    lines(survfit(Surv(out2$Times[, k], out2$event[, k]) ~ 1), col = "lightgrey",
          conf.int = FALSE, fun = "event")
}
lines(survfit(Surv(Times, event) ~ 1), fun = "event")
legend("bottomright", c("replicated data", "observed data"), lty = 1,
       col = c("lightgrey", "black"), bty = "n", cex = 0.8)


xvals <- seq(0, max(Times), length.out = 500)
plot(range(Times), c(0, 1), type = "n", xlab = "Time",
     ylab = "Cumulative Incidence")
for (k in seq_len(ncol(out3$Times))) {
    lines(survfit(Surv(out3$Times[, k], out3$event[, k]) ~ 1), col = "lightgrey",
          conf.int = FALSE, fun = "event")
}
lines(survfit(Surv(Times, event) ~ 1), fun = "event")
legend("bottomright", c("replicated data", "observed data"), lty = 1,
       col = c("lightgrey", "black"), bty = "n", cex = 0.8)

simulate1.jm <-
    function (object, nsim = 1L, seed = NULL,
              process = c("longitudinal", "event"),
              random_effects = c("posterior_means", "mcmc", "prior"),
              Fforms_fun = NULL, tol = 0.001, ...) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1L)
        if (is.null(seed))
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
        process <- match.arg(process)
        random_effects <- match.arg(random_effects)
        # information from fitted joint model
        n <- object$model_data$n
        idL_lp <- object$model_data$idL_lp
        ind_RE <- object$model_data$ind_RE
        X <- object$model_data$X
        Z <- object$model_data$Z
        n_outcomes <- length(idL_lp)
        families <- object$model_info$families
        has_sigmas <- as.logical(object$model_data$has_sigmas)
        Bspline_dgr <- object$control$Bsplines_degree
        knots <- object$control$knots[[1]]
        W <- object$model_data$W_h
        Times <- object$model_data$Time_right
        event <- object$model_data$delta
        dataS <- object$model_data$dataS

        # MCMC results
        ncz <- sum(sapply(Z, ncol))
        ind_betas <- grep("betas", names(object$statistics$Mean), fixed = TRUE)
        mcmc_betas <- object$mcmc[ind_betas]
        mcmc_betas[] <- lapply(mcmc_betas, function (x) do.call('rbind', x))
        mcmc_sigmas <- matrix(0.0, nrow(mcmc_betas[[1]]), n_outcomes)
        mcmc_sigmas[, has_sigmas] <- do.call('rbind', object$mcmc$sigmas)
        mcmc_bs_gammas <- do.call('rbind', object$mcmc$bs_gammas)
        has_gammas <- !is.null(object$mcmc$gammas)
        if (has_gammas) mcmc_gammas <- do.call('rbind', object$mcmc$gammas)
        mcmc_W_std_gammas <- do.call('rbind', object$mcmc$W_std_gammas)
        mcmc_alphas <- do.call('rbind', object$mcmc$alphas)
        mcmc_Wlong_std_alphas <- do.call('rbind', object$mcmc$Wlong_std_alphas)
        # random effects
        b <- ranef(object)
        if (random_effects == "mcmc") {
            mcmc_RE <- dim(object$mcmc[["b"]][[1L]])[3L] > 1L
            if (mcmc_RE) {
                mcmc_b <- abind::abind(object$mcmc[["b"]])
            } else {
                stop("refit the model using 'jm(..., save_random_effects = TRUE)'.\n")
            }
        }
        get_D <- function (x, n) {
            m <- matrix(0.0, n, n)
            m[lower.tri(m, TRUE)] <- x
            m <- m + t(m)
            diag(m) <- diag(m) * 0.5
            m
        }
        xx <- do.call('rbind', object$mcmc$D)
        mcmc_D <- array(0.0, c(ncz, ncz, nrow(xx)))
        for (m in seq_len(nrow(mcmc_betas[[1L]]))) {
            mcmc_D[, , m] <- get_D(xx[m, ], ncz)
        }
        # simulate outcome vectors
        if (process == "longitudinal") {
            sim_fun <- function (family, n, mu, phi) {
                switch(family$family,
                       "gaussian" = rnorm(n, mu, phi),
                       "Student's-t" = mu + phi * rt(n, df = 4),
                       "binomial" = rbinom(n, 1, mu),
                       "poisson" = rpois(n, mu),
                       "beta" = rbeta(n, shape1 = mu * phi, shape2 = phi * (1.0 - mu)))
            }
            val <- vector("list", nsim)
            indices <- sample(nrow(mcmc_betas[[1]]), nsim)
            for (j in seq_len(nsim)) {
                # parameters
                jj <- indices[j]
                betas <- lapply(mcmc_betas, function (x) x[jj, ])
                sigmas <- mcmc_sigmas[jj, ]
                bb <-
                    switch(random_effects,
                           "mcmc" = mcmc_b[, , jj],
                           "posterior_means" = b,
                           "prior" = MASS::mvrnorm(n, rep(0, ncz), mcmc_D[, , jj]))
                y <- vector("list", n_outcomes)
                for (i in seq_len(n_outcomes)) {
                    FE <- c(X[[i]] %*% betas[[i]])
                    RE <- rowSums(Z[[i]] * bb[idL_lp[[i]], ind_RE[[i]]])
                    eta <- FE + RE
                    mu <- families[[i]]$linkinv(eta)
                    y[[i]] <- sim_fun(families[[i]], length(mu), mu, sigmas[i])
                }
                val[[j]] <- y
            }
            names(val) <- paste0("sim_", seq_len(nsim))
        } else {
            if (is.null(Fforms_fun) || !is.function(Fforms_fun)) {
                stop("you need to provide the 'Fforms_fun' function; ",
                     "see the examples in `?simulate.jm` for more information.\n")
            }
            invS <- function (t, log_u, subj, rescale_factor) {
                hazard <- function (t, subj, rescale_factor) {
                    W0 <- splineDesign(knots, t, Bspline_dgr + 1L, outer.ok = TRUE)
                    log_h0 <- c(W0 %*% bs_gammas) - rescale_factor
                    covariates <- if (has_gammas) c(W[subj, , drop = FALSE] %*% gammas) else 0.0
                    long <- c(Fforms_fun(t, betas, b[subj, ], dataS[subj, ]) %*% alphas)
                    exp(log_h0 + covariates + long)
                }
                integrate(hazard, lower = 0, upper = t, subj = subj,
                          rescale_factor = rescale_factor, rel.tol = tol)$value + log_u
            }
            valT <- eventT <- matrix(0.0, n, nsim)
            indices <- sample(nrow(mcmc_betas[[1]]), nsim)
            for (j in seq_len(nsim)) {
                jj <- indices[j]
                betas <- lapply(mcmc_betas, function (x) x[jj, ])
                bb <-
                    switch(random_effects,
                           "mcmc" = mcmc_b[, , jj],
                           "posterior_means" = b,
                           "prior" = MASS::mvrnorm(n, rep(0, nRE), mcmc_D[, , jj]))
                bs_gammas <- mcmc_bs_gammas[jj, ]
                if (has_gammas) gammas <- mcmc_gammas[jj, ]
                alphas <- mcmc_alphas[jj, ]
                Wlong_std_alphas <- mcmc_Wlong_std_alphas[jj, ]
                W_std_gammas <- mcmc_W_std_gammas[jj, ]
                rescale_factor <- Wlong_std_alphas + W_std_gammas
                Up <- max(Times) * 1.05
                rep_Times <- numeric(n)
                if (j == 1) {
                    invS(0.5, runif(1L), 1, rescale_factor)
                }
                for (i in seq_len(n)) {
                    log_u <- log(runif(1L))
                    Root <-
                        try(uniroot(invS, interval = c(0, Up), log_u = log_u,
                                    subj = i, rescale_factor = rescale_factor,
                                    tol = tol, f.lower = log_u)$root,
                            silent = TRUE)
                    rep_Times[i] <- if (!inherits(Root, "try-error")) Root else Up
                }
                valT[, j] <- rep_Times
                eventT[, j] <- as.numeric(rep_Times < Up)
            }
            colnames(valT) <- colnames(eventT) <- paste0("sim_", seq_len(nsim))
            val <- list(Times = valT, event = eventT)
        }
        attr(val, "seed") <- RNGstate
        val
    }


simulate2.jm <-
    function (object, nsim = 1L, seed = NULL,
              process = c("longitudinal", "event"),
              random_effects = c("posterior_means", "mcmc", "prior"),
              Fforms_fun = NULL, tol = 0.001, subdivisions = 2L, ...) {
        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            runif(1L)
        if (is.null(seed))
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }
        process <- match.arg(process)
        random_effects <- match.arg(random_effects)
        # information from fitted joint model
        n <- object$model_data$n
        idL_lp <- object$model_data$idL_lp
        ind_RE <- object$model_data$ind_RE
        X <- object$model_data$X
        Z <- object$model_data$Z
        n_outcomes <- length(idL_lp)
        families <- object$model_info$families
        has_sigmas <- as.logical(object$model_data$has_sigmas)
        Bspline_dgr <- object$control$Bsplines_degree
        knots <- object$control$knots[[1]]
        W <- object$model_data$W_h
        Times <- object$model_data$Time_right
        event <- object$model_data$delta
        dataS <- object$model_data$dataS

        # MCMC results
        ncz <- sum(sapply(Z, ncol))
        ind_betas <- grep("betas", names(object$statistics$Mean), fixed = TRUE)
        mcmc_betas <- object$mcmc[ind_betas]
        mcmc_betas[] <- lapply(mcmc_betas, function (x) do.call('rbind', x))
        mcmc_sigmas <- matrix(0.0, nrow(mcmc_betas[[1]]), n_outcomes)
        if (has_sigmas) mcmc_sigmas[, has_sigmas] <- do.call('rbind', object$mcmc$sigmas)
        mcmc_bs_gammas <- do.call('rbind', object$mcmc$bs_gammas)
        has_gammas <- !is.null(object$mcmc$gammas)
        if (has_gammas) mcmc_gammas <- do.call('rbind', object$mcmc$gammas)
        mcmc_W_std_gammas <- do.call('rbind', object$mcmc$W_std_gammas)
        mcmc_alphas <- do.call('rbind', object$mcmc$alphas)
        mcmc_Wlong_std_alphas <- do.call('rbind', object$mcmc$Wlong_std_alphas)
        # random effects
        b <- ranef(object)
        if (random_effects == "mcmc") {
            mcmc_RE <- dim(object$mcmc[["b"]][[1L]])[3L] > 1L
            if (mcmc_RE) {
                mcmc_b <- abind::abind(object$mcmc[["b"]])
            } else {
                stop("refit the model using 'jm(..., save_random_effects = TRUE)'.\n")
            }
        }
        get_D <- function (x, n) {
            m <- matrix(0.0, n, n)
            m[lower.tri(m, TRUE)] <- x
            m <- m + t(m)
            diag(m) <- diag(m) * 0.5
            m
        }
        xx <- do.call('rbind', object$mcmc$D)
        mcmc_D <- array(0.0, c(ncz, ncz, nrow(xx)))
        for (m in seq_len(nrow(mcmc_betas[[1L]]))) {
            mcmc_D[, , m] <- get_D(xx[m, ], ncz)
        }
        # simulate outcome vectors
        if (process == "longitudinal") {
            sim_fun <- function (family, n, mu, phi) {
                switch(family$family,
                       "gaussian" = rnorm(n, mu, phi),
                       "Student's-t" = mu + phi * rt(n, df = 4),
                       "binomial" = rbinom(n, 1, mu),
                       "poisson" = rpois(n, mu),
                       "negative binomial" = rnbinom(n, size = phi, mu = mu),
                       "beta" = rbeta(n, shape1 = mu * phi, shape2 = phi * (1.0 - mu)))
            }
            val <- vector("list", nsim)
            indices <- sample(nrow(mcmc_betas[[1]]), nsim)
            for (j in seq_len(nsim)) {
                # parameters
                jj <- indices[j]
                betas <- lapply(mcmc_betas, function (x) x[jj, ])
                sigmas <- mcmc_sigmas[jj, ]
                bb <-
                    switch(random_effects,
                           "mcmc" = mcmc_b[, , jj],
                           "posterior_means" = b,
                           "prior" = MASS::mvrnorm(n, rep(0, ncz), mcmc_D[, , jj]))
                y <- vector("list", n_outcomes)
                for (i in seq_len(n_outcomes)) {
                    FE <- c(X[[i]] %*% betas[[i]])
                    RE <- rowSums(Z[[i]] * bb[idL_lp[[i]], ind_RE[[i]]])
                    eta <- FE + RE
                    mu <- families[[i]]$linkinv(eta)
                    y[[i]] <- sim_fun(families[[i]], length(mu), mu, sigmas[i])
                }
                val[[j]] <- y
            }
            names(val) <- paste0("sim_", seq_len(nsim))
        } else {
            if (is.null(Fforms_fun) || !is.function(Fforms_fun)) {
                stop("you need to provide the 'Fforms_fun' function; ",
                     "see the examples in `?simulate.jm` for more information.\n")
            }
            if (length(unique(object$model_data$strata)) > 1L) {
                stop("'simulate.jm()' does not currently support stratified joint models.\n")
            }
            hazard <- function (time, log_u, subj, rescale_factor) {
                W0 <- splineDesign(knots, time, Bspline_dgr + 1L, outer.ok = TRUE)
                log_h0 <- c(W0 %*% bs_gammas) - rescale_factor
                covariates <- if (has_gammas) c(W[subj, , drop = FALSE] %*% gammas) else 0.0
                long <- c(Fforms_fun(time, betas, b[subj, ], dataS[subj, ]) %*% alphas)
                exp(log_h0 + covariates + long)
            }
            invS <- function (time, log_u, subj, rescale_factor) {
                integrate(hazard, lower = 0.0, upper = time, subj = subj,
                          rescale_factor = rescale_factor, rel.tol = tol,
                          subdivisions = subdivisions, stop.on.error = FALSE)$value + log_u
            }
            nr_root <- function (interval, fn, gr, ..., tol = tol, iter = 35L) {
                Low <- low <- interval[1L]
                Upp <- upp <- interval[2L]
                fn_Low <- fn(Low, ...)
                fn_Upp <- fn(Upp, ...)
                if (fn_Upp < 0) return(Upp)
                tt <- tt_old <- sum(interval) / 2
                for (i in seq_len(iter)) {
                    ffn <- fn(tt, ...)
                    # check convergence
                    if (abs(ffn) < tol) return(tt)
                    # if proposed value outside interval do bisection,
                    # else Newton-Raphson
                    ggr <- gr(tt, ...)
                    tt <- tt - ffn / ggr
                    if (tt < Low || tt > Upp) {
                        if (ffn < 0 && ffn > fn_Low) {
                            low <- tt_old
                            fn_Low <- ffn
                        } else {
                            upp <- tt_old
                        }
                        tt <- tt_old <- (low + upp) / 2
                    }
                }
                tt
            }
            valT <- eventT <- matrix(0.0, n, nsim)
            indices <- sample(nrow(mcmc_betas[[1]]), nsim)
            for (j in seq_len(nsim)) {
                jj <- indices[j]
                betas <- lapply(mcmc_betas, function (x) x[jj, ])
                bb <-
                    switch(random_effects,
                           "mcmc" = mcmc_b[, , jj],
                           "posterior_means" = b,
                           "prior" = MASS::mvrnorm(n, rep(0, nRE), mcmc_D[, , jj]))
                bs_gammas <- mcmc_bs_gammas[jj, ]
                if (has_gammas) gammas <- mcmc_gammas[jj, ]
                alphas <- mcmc_alphas[jj, ]
                Wlong_std_alphas <- mcmc_Wlong_std_alphas[jj, ]
                W_std_gammas <- mcmc_W_std_gammas[jj, ]
                rescale_factor <- Wlong_std_alphas + W_std_gammas
                Up <- max(Times) * 1.05
                rep_Times <- numeric(n)
                if (j == 1) {
                    invS(0.5, -0.69, 1, rescale_factor)
                }
                for (i in seq_len(n)) {
                    log_u <- log(runif(1L))
                    rep_Times[i] <-
                        nr_root(c(0, Up), invS, hazard, subj = i, log_u = log_u,
                                rescale_factor = rescale_factor, tol = tol)
                }
                valT[, j] <- rep_Times
                eventT[, j] <- as.numeric(rep_Times < Up)
            }
            colnames(valT) <- colnames(eventT) <- paste0("sim_", seq_len(nsim))
            val <- list(Times = valT, event = eventT)
        }
        attr(val, "seed") <- RNGstate
        val
    }


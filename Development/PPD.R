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
fm2 <- mixed_model(ascites ~ year, data = pbc2, random = ~ year | id,
                   family = binomial())

# the joint model
jointFit1 <- jm(CoxFit, list(fm1), time_var = "year", save_random_effects = TRUE)


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


ND <- pbc2[pbc2$id == 2, ]
prs <- predict(jointFit, newdata = ND, return_params_mcmc = TRUE)
str(prs$params_mcmc)

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

# Cox model for the composite event death or transplantation
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
pbc2$status2 <- as.numeric(pbc2$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)

# a linear mixed model for log serum bilirubin
fm1 <- lme(log(serBilir) ~ ns(year, 3) * sex, data = pbc2,
           random = list(id = pdDiag(~ ns(year, 3))))
fm2 <- lme(log(serBilir) ~ poly(year, 2) * sex, data = pbc2,
           random = list(id = pdDiag(~ poly(year, 2))))


# the joint model
jointFit1 <- jm(CoxFit, list(fm1), time_var = "year", save_random_effects = TRUE)
jointFit2 <- jm(CoxFit, list(fm2), time_var = "year", save_random_effects = TRUE)


MISE_perid <- function (obs, reps, id, times, t0, Dt) {
    trapezoid_rule <- function (f, x) {
        sum(0.5 * diff(x) * (f[-length(x)] + f[-1L]))
    }
    loess.smooth2 <- function (x, y) {
        loess.smooth(x, y, degree = 2, span = 0.75,
                     family = "gaussian", evaluation = 200)
    }
    smooth <- function (x, y) {
        n <- length(x)
        if (n > 5) {
            loess.smooth2(x, y)
        } else if (n > 1 && n <= 5) {
            approx(x, y)
        } else {
            list(x = NA_real_, y = NA_real_)
        }
    }
    gof_fun <- function (y, times, id, type) {
        if (type == "variogram") {
            ls <- loess.smooth2(times, y)
            ind <- findInterval(times, ls$x)
            rr <- y - ls$y[ind]
            variogram(rr, times, id)[[1L]]
        } else if (type == "variance-function") {
            ls <- loess.smooth2(times, y)
            ind <- findInterval(times, ls$x)
            rr <- y - ls$y[ind]
            sigma <- sqrt(sum(rr * rr) / (length(rr) - 5.35))
            rr <- sqrt(abs(rr / sigma))
            cbind(times, rr)
        } else {
            cbind(times, y)
        }
    }
    Obs_ave <- gof_fun(obs, times, id, "average")
    Obs_vario <- gof_fun(obs, times, id, "variogram")
    F_obs_ave <- mapply(smooth, y = split(Obs_ave[, 2L], id),
                        x = split(Obs_ave[, 1L], id), SIMPLIFY = FALSE)
    ni <- tapply(id, id, length)
    id_long <- rep(unique(id), sapply(ni, function (n) ncol(combn(n, 2))))
    F_obs_vario <- mapply(smooth, y = split(Obs_vario[, 2L], id_long),
                          x = split(Obs_vario[, 1L], id_long), SIMPLIFY = FALSE)
    mise <- function (obs, rep, t0, Dt, type) {
        xx <- obs$x
        ind <- if (type == "average") {
            xx > t0 - Dt & xx < t0
        } else {
            xx > 0 & xx < Dt
        }
        trapezoid_rule((obs$y[ind] - rep$y[ind])^2, xx[ind])
    }
    n <- length(unique(id))
    M <- ncol(reps)
    MISE <- matrix(0.0, n, M)
    for (m in seq_len(M)) {
        reps_ave <- gof_fun(reps[, m], times, id, "average")
        reps_vario <- gof_fun(reps[, m], times, id, "variogram")
        F_reps_ave <-
            mapply(smooth, y = split(reps_ave[, 2L], id),
                   x = split(reps_ave[, 1L], id), SIMPLIFY = FALSE)
        F_reps_vario <-
            mapply(smooth, y = split(reps_vario[, 2L], id_long),
                   x = split(reps_vario[, 1L], id_long), SIMPLIFY = FALSE)
        MISE[, m] <- mapply(mise, obs = F_obs_ave, rep = F_reps_ave,
                            MoreArgs = list(t0 = t0, Dt = Dt, type = "average")) +
            200 * mapply(mise, obs = F_obs_vario, rep = F_reps_vario,
                   MoreArgs = list(t0 = t0, Dt = Dt, type = "variogram"))
    }
    rowMeans(MISE)
}

t0 <- 6
Dt <- 2
ND <- pbc2[ave(pbc2$year, pbc2$id, FUN = max) > t0, ]
ND$id <- match(ND$id, unique(ND$id))
ND_before <- ND[ND$year <= t0, ]
ND_after <- ND[ND$year > t0, ]

#obs = sims1$outcome[[1]]
#reps = sims1[[1L]]
#id = ND_before$id
#times = ND_before$year

prs1 <- predict(jointFit1, newdata = ND_before, newdata2 = ND_after,
                return_params_mcmc = TRUE)
sims1 <- simulate(jointFit1, nsim = 200L, newdata = ND_before,
                  include_outcome = TRUE, random_effects = "mcmc",
                  params_mcmc = prs1$newdata$params_mcmc)
MISEs1 <- MISE_perid(sims1$outcome[[1]], sims1[[1L]], ND_before$id,
                     ND_before$year, t0, Dt)
#
prs2 <- predict(jointFit2, newdata = ND_before, newdata2 = ND_after,
                return_params_mcmc = TRUE)
sims2 <- simulate(jointFit2, nsim = 200L, newdata = ND_before,
                  include_outcome = TRUE, random_effects = "mcmc",
                  params_mcmc = prs2$newdata$params_mcmc)
MISEs2 <- MISE_perid(sims2$outcome[[1]], sims2[[1L]], ND_before$id,
                     ND_before$year, t0, Dt)


ppcheck(jointFit1, nsim = 500L, newdata = ND_before[ND_before$id == 4, ],
        random_effects = "mcmc", params_mcmc = prs1$params_mcmc,
        type = "average", CI_loess = TRUE)
ppcheck(jointFit2, nsim = 500L, newdata = ND_before[ND_before$id == 4, ],
        random_effects = "mcmc", params_mcmc = prs2$params_mcmc,
        type = "average", CI_loess = TRUE)



out <- data.frame(
    id = unique(ND_before$id),
    MISEs1 = MISEs1, MISEs2 = MISEs2,
    Model = ifelse(MISEs1 < MISEs2, "Model I", "Model II"))

Preds <- data.frame(id = ND_after$id, Obs = log(ND_after$serBilir),
                    Preds1 = prs1$newdata2$preds$serBilir,
                    Preds2 = prs2$newdata2$preds$serBilir)
Preds$Opt_Model <- out$Model[ND_after$id]
Preds$Opt_Preds <- with(Preds, ifelse(Opt_Model == "Model I", Preds1, Preds2))
Preds$MISE1 <- unlist(mapply(function (obs, pred) rep(mean((obs - pred)^2), length(obs)),
                      obs = split(Preds$Obs, Preds$id),
                      pred = split(Preds$Preds1, Preds$id)))
Preds$MISE2 <- unlist(mapply(function (obs, pred) rep(mean((obs - pred)^2), length(obs)),
                             obs = split(Preds$Obs, Preds$id),
                             pred = split(Preds$Preds2, Preds$id)))
Preds$MISE_opt <- unlist(mapply(function (obs, pred) rep(mean((obs - pred)^2), length(obs)),
                                obs = split(Preds$Obs, Preds$id),
                                pred = split(Preds$Opt_Preds, Preds$id)))
Preds$Selection <- with(Preds, ifelse(MISE_opt == pmin(MISE1, MISE2), "correct", "wrong"))

with(Preds, mean((Obs - Preds1)^2))
with(Preds, mean((Obs - Preds2)^2))
with(Preds, mean((Obs - Opt_Preds)^2))


##########################################################################################
# Author: D. Rizopoulos                                                                  #
# Aim: Test function jm()                                                                #
#########################################################################################
library("survival")
library("nlme")
library("GLMMadaptive")
library("coda")
library("lattice")
library("splines")
library("matrixStats")
source("./R/jm.R")
source("./R/jm_fit.R")
source("./R/help_functions.R")
source("./R/basic_methods.R")
source("./R/create_Wlong_mats.R")
source("./Development/jm/PBC_data.R")
Rcpp::sourceCpp('src/mcmc_fit.cpp')

simulateJoint <- function (alpha = 0.0, Dalpha = 0, n = 500,
                           mean.Cens = 35) {
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
    trueValues <- list(betas = betas, tau = 1/sigma.y^2, gammas = gammas,
                       alphas = alpha, Dalphas = Dalpha, sigma.t = phi,
                       inv.D = solve(D), b = b)

    # return list
    list(DF = dat, DF.id = dat.id, trueValues = trueValues)
}

Data <- simulateJoint()
lmeFit <- lme(y ~ ns(time, 3), data = Data$DF,
              random = list(id = pdDiag(form = ~ ns(time, 3))),
              control = lmeControl(opt = "optim", niterEM = 45))
coxFit <- coxph(Surv(Time, event) ~ group + age, data = Data$DF.id)

obj <- jm(coxFit, list(lmeFit), time_var = "time", priors = list(sigmas_shape = 0.1))

summary(obj)
#coda::traceplot(obj$mcmc$D)
coda::autocorr.diag(obj$mcmc$D)
coda::cumuplot(obj$mcmc$alphas)

traceplot(obj)
ggtraceplot(obj)
gelman_diag(obj)

#obj
model_data <- obj$model_data
model_info <- obj$model_info
initial_values <- obj$initial_values
priors <- obj$priors
vcov_prop <- obj$vcov_prop
control <- obj$control

Surv_object = coxFit
Mixed_objects = lmeFit
time_var = "time"
functional_forms = NULL
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con
control$n_chains = 1


#
pbc2$prothrombin[pbc2$id == levels(pbc2$id)[1L]] <- NA
pbc2$prothrombin[pbc2$id == levels(pbc2$id)[2L]] <- NA

fm1 <- lme(log(serBilir) ~ ns(year, 2, B = c(0, 9)) * sex, data = pbc2,
           random = ~ ns(year, 2, B = c(0, 9)) | id)

fm2 <- lme(prothrombin ~ ns(year, 2, B = c(0, 9)) * sex, data = pbc2,
           random = ~ ns(year, 2, B = c(0, 9)) | id,
           na.action = na.exclude, control = lmeControl(opt = "optim"))

fm3 <- mixed_model(ascites ~ year * sex, data = pbc2, random = ~ year | id,
                   family = binomial())
Mixed <- list(fm1, fm2, fm3)
Cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

system.time(obj <- jm(Cox, Mixed, time_var = "year"))

summary(obj)
traceplot(obj, "alphas")
gelman_diag(obj)

Surv_object = Cox
Mixed_objects = Mixed
time_var = 'year'
functional_forms = NULL
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con
control$n_chains = 1



##########################################################################################
##########################################################################################

fm1 <- lme(log(serBilir) ~ year * sex + I(year^2) + age + prothrombin,
           data = pbc2, random = ~ year | id)
fm2 <- lme(serChol ~ ns(year, 3) + sex + age, data = pbc2, random = ~ year | id,
           na.action = na.exclude)
fm3 <- mixed_model(hepatomegaly ~ sex + age + year, data = pbc2,
                   random = ~ year | id, family = binomial())
fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ year | id, family = binomial())

CoxFit <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

#CoxFit <- survreg(Surv(years, yearsU, status3, type = "interval") ~ 1,
#                  data = pbc2.id, model = TRUE)

fForms <- list("log(serBilir)" = ~ value(log(serBilir)) + slope(log(serBilir)) +
                   value(log(serBilir)):sex,
               "serChol" = ~ value(serChol) + slope(serChol),
               "hepatomegaly" = ~ value(hepatomegaly) + sex,
               "ascites" = ~ value(ascites) + area(ascites))

test <- jm(CoxFit, list(fm1, fm2, fm3, fm4), time_var = "year",
           functional_forms = fForms)

####

Surv_object = CoxFit
Mixed_objects = list(fm1, fm2, fm3, fm4)
time_var = "year"
functional_forms = NULL
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL


################################################################################
################################################################################


fm <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
CoxFit <- coxph(Surv(years, status2) ~ age, data = pbc2.id)
ff <- list("log(serBilir)" = ~ area(log(serBilir)))
test <- jm(CoxFit, fm, time_var = "year", functional_forms = ff)
summary(test)


####

Surv_object = CoxFit
Mixed_objects = fm
time_var = "year"
functional_forms = NULL
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con
control$n_chains = 1

################################################################################

fm <- lme(CD4 ~ obstime * drug, data = aids, random = ~ obstime | patient)
gm <- coxph(Surv(Time, death) ~ drug, data = aids.id)
jmFit <- jm(gm, fm, time_var = "obstime")
summary(jmFit)











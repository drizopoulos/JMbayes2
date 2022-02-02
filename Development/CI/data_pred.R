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
source("./R/predict_funs.R")
load("./Development/CI/model.RData")
source("./Development/CI/prepare_data.R")
Rcpp::sourceCpp('src/mcmc_fit.cpp')


newdataL <- pbc2[pbc2$id %in% c(2, 81), ]
newdataE <- pbc2_CR[pbc2_CR$id %in% c(2, 81), ]
newdataE$event <- 0
newdataE$stop <- c(9, 9, 3.175994, 7, 3.175994, 7)
newdata <- list(newdataL = newdataL, newdataE = newdataE)
##
newdataL2 <- tail(newdataL, 1)
newdataL2$IE <- 1
newdataL2$S <- newdataL2$year
newdataE2 <- newdataE
newdataE2$IE <- 1
newdata2 <- list(newdataL = newdataL2, newdataE = newdataE2)

ND_long <- pbc2[pbc2$id %in% c(12, 81), ]
ND_event <- pbc2.idCR[pbc2.idCR$id %in% c(12, 81), ]
ND_event$status2 <- 0
ND <- list(newdataL = ND_long, newdataE = ND_event)

object = jointFit3
newdata = ND
newdata2 = NULL #newdata2
times = c(7.1, 8.1, 9.1)
process = "event"
type_pred = "response"
type = "subject_specific"
level = 0.95; return_newdata = TRUE
n_samples = 200L; n_mcmc = 55L; cores = NULL
seed = 123L

predLong2 <- predict(jointFit, newdata = newdata,
                     times = seq(5, 12, length.out = 51),
                     return_newdata = TRUE)
predSurv <- predict(jointFit, newdata = newdata, process = "event",
                    return_newdata = TRUE)

plot(predSurv)

plot(predLong2, predSurv)



##############

newdata <- pbc2[pbc2$id %in% c(14, 2), ]
newdata$years[newdata$id == 2] <- 3.3
newdata$years[newdata$id == 14] <- 9
newdata$status2 <- 0

object = jointFit2
newdata = newdata
newdata2 = NULL
times = NULL
process = "event"
type_pred = "response"
type = "subject_specific"
level = 0.95; return_newdata = TRUE; return_mcmc = FALSE
n_samples = 200L; n_mcmc = 55L; cores = NULL
seed = 123L

##############


t0 <- 5
ND <- pbc2[pbc2$id %in% c(25, 93), ]
ND <- ND[ND$year < t0, ]
ND$status2 <- 0
ND$years <- t0
newdata <- ND

predLong2 <- predict(jointFit2, newdata = ND,
                     times = seq(t0, 12, length.out = 51),
                     return_newdata = TRUE)
predSurv <- predict(jointFit2, newdata = ND, process = "event",
                    times = seq(t0, 12, length.out = 51),
                    return_newdata = TRUE)

plot(predLong2, predSurv)



pbc2.id$event <- as.numeric(pbc2.id$status != "alive")
CoxFit <- coxph(Surv(years, event) ~ age, data = pbc2.id)
fm1 <- lme(log(serBilir) ~ ns(year, 3) * sex, data = pbc2,
           random = ~ ns(year, 3) | id, control = lmeControl(opt = 'optim'))

fm2 <- lme(prothrombin ~ ns(year, 2) * sex, data = pbc2,
           random = ~ ns(year, 2) | id, control = lmeControl(opt = 'optim'))

fm3 <- mixed_model(ascites ~ year * sex, data = pbc2,
                   random = ~ year | id, family = binomial())

jointFit <- jm(CoxFit, list(fm1, fm2, fm3), time_var = "year")


prs <- predict(jointFit, newdata, process = "event")


fm1 <- lme(log(serBilir) ~ poly(year, 2) * drug, data = pbc2,
           random = ~ poly(year, 2) | id)

CR_forms <- list(
    "log(serBilir)" = ~ value(log(serBilir)):CR
)
CoxFit_CR <- coxph(Surv(years, status2) ~ (age + drug) * strata(CR),
                   data = pbc2.idCR)

jointFit4 <- jm(CoxFit_CR, list(fm1), time_var = "year",
                functional_forms = CR_forms)


ND_long <- pbc2[pbc2$id == 81, ]
ND_event <- pbc2.idCR[pbc2.idCR$id == 81, ]
ND_event$status2 <- 0
ND <- list(newdataL = ND_long, newdataE = ND_event)

predLong <- predict(jointFit3, newdata = ND, return_newdata = TRUE,
                    times = seq(6.5, 15, length = 25))

predEvent <- predict(jointFit3, newdata = ND, return_newdata = TRUE,
                     process = "event")

plot(predLong, predEvent, outcomes = 1:2, pos_ylab_long = c(0.1, 20))
plot(predLong, predEvent, outcomes = 1, pos_ylab_long = c(0.1, 20))


predEvent <- predict(jointFit3, newdata = ND, process = "event",
                     return_mcmc = T)

rowSums(matrix(predEvent$mcmc[, 21], nc = 2))


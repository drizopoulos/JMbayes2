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


newdataL <- pbc2[pbc2$id %in% c(14, 13), ]
newdataE <- pbc2_CR[pbc2_CR$id %in% c(14, 13), ]
newdataE$event <- 0
newdata <- list(newdataL = newdataL, newdataE = newdataE)

object = jointFit
newdata = newdata
newdata2 = NULL
times = NULL
process = "event"
type_pred = "response"
type = "subject_specific"
level = 0.95; return_newdata = TRUE
n_samples = 200L; n_mcmc = 55L; cores = NULL
seed = 123L

predLong2 <- predict(jointFit, newdata = ND,
                     times = seq(t0, 12, length.out = 51),
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
level = 0.95; return_newdata = TRUE
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



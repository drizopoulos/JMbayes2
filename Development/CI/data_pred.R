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

ii <- 1:5
ND_long <- pbc2[pbc2$id == 25, ][ii, ]
ND_event <- pbc2.idCR[pbc2.idCR$id == 25, ]
ND_event$status2 <- 0
ND_event$years <- max(ND_long[["year"]])
ND <- list(newdataL = ND_long, newdataE = ND_event)

predLong <- predict(jointFit4, newdata = ND, return_newdata = TRUE,
                    times = seq(0.5, 25, length = 25))

predEvent <- predict(jointFit4, newdata = ND, return_newdata = TRUE,
                     process = "event")

#plot(predLong, predEvent, outcomes = 1:2, pos_ylab_long = c(0.1, 20))
plot(predLong, predEvent, outcomes = 1, pos_ylab_long = c(0.1, 20))


predEvent <- predict(jointFit3, newdata = ND, process = "event",
                     return_mcmc = T)

rowSums(matrix(predEvent$mcmc[, 21], nc = 2))


###############################################################################

fm1 <- lme(log(serBilir) ~ ns(year, 3) * sex, data = pbc2,
           random = ~ ns(year, 3) | id, control = lmeControl(opt = 'optim'))

fm2 <- lme(prothrombin ~ ns(year, 2) * sex, data = pbc2,
           random = ~ ns(year, 2) | id, control = lmeControl(opt = 'optim'))

fm3 <- mixed_model(ascites ~ year * sex, data = pbc2,
                   random = ~ year | id, family = binomial())
pbc2.id$event <- as.numeric(pbc2.id$status != "alive")
CoxFit <- coxph(Surv(years, event) ~ drug + age, data = pbc2.id)
jointFit <- jm(CoxFit, list(fm1, fm2, fm3), time_var = "year")

t0 <- 5
ND <- pbc2[pbc2$id %in% c(25, 93), ]
ND <- ND[ND$year < t0, ]
ND$status2 <- 0
ND$years <- t0


object = jointFit
newdata = ND
newdata2 = NULL
times = NULL
process = "event"
type_pred = "response"
type = "subject_specific"
level = 0.95; return_newdata = TRUE; return_mcmc = FALSE
n_samples = 200L; n_mcmc = 55L; cores = NULL
seed = 123L


predLong1 <- predict(jointFit, newdata = ND, return_newdata = TRUE)
plot(predLong1)
predLong2 <- predict(jointFit, newdata = ND,
                     times = seq(t0, 12, length.out = 51),
                     return_newdata = TRUE)
plot(predLong2, outcomes = 2, subject = 93)

predSurv <- predict(jointFit, newdata = ND, process = "event",
                    times = seq(t0, 12, length.out = 51),
                    return_newdata = TRUE)

plot(predLong2, predSurv)

cols <- c('#F25C78', '#D973B5', '#F28322')
plot(predLong2, predSurv, outcomes = 1:3, subject = 93,
     fun_long = list(exp, identity, identity),
     fun_event = function (x) 1 - x,
     ylab_event = "Survival Probabilities",
     ylab_long = c("Serum Bilirubin", "Prothrombin", "Ascites"),
     bg = '#132743', col_points = cols, col_line_long = cols,
     col_line_event = '#F7F7FF', col_axis = "white",
     fill_CI_long = c("#F25C7880", "#D973B580", "#F2832280"),
     fill_CI_event = "#F7F7FF80",
     pos_ylab_long = c(1.9, 1.9, 0.08))


x = predLong2; x2 = predSurv; subject = 1; outcomes = 1;
fun_long = NULL; fun_event = NULL
CI_long = TRUE; CI_event = TRUE;
xlab = "Follow-up Time"; ylab_long = NULL
ylab_event = "Cumulative Risk"; main = "";
lwd_long = 2; lwd_event = 2
ylim_long_outcome_range = TRUE
col_line_long = "#0000FF"
col_line_event = c("#FF0000", "#03BF3D", "#8000FF")
pch_points = 16; col_points = "blue"; cex_points = 1;
fill_CI_long = "#0000FF4D"
fill_CI_event = c("#FF00004D", "#03BF3D4D", "#8000FF4D")
cex_xlab = 1; cex_ylab_long = 1; cex_ylab_event = 1
cex_main = 1; cex_axis = 1; col_axis = "black"
pos_ylab_long = c(0.1, 2, 0.08); bg = "white"





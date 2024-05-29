###############################################################################
################################################################################
################################################################################

if (FALSE) {
    library("JMbayes2")
    pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
    CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
    fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
    #fm2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)
    #fm3 <- mixed_model(ascites ~ year + sex, data = pbc2,
     #                  random = ~ year | id, family = binomial())

    jointFit1 <- jm(CoxFit, list(fm1), time_var = "year")

    pbc2.idCR <- crisk_setup(pbc2.id, statusVar = "status", censLevel = "alive",
                             nameStrata = "CR")
    CoxFit_CR <- coxph(Surv(years, status2) ~ (age + drug) * strata(CR),
                       data = pbc2.idCR)
    fm1 <- lme(log(serBilir) ~ poly(year, 2) * drug, data = pbc2,
               random = ~ poly(year, 2) | id)
    fm2 <- lme(prothrombin ~ year * drug, data = pbc2, random = ~ year | id)
    CR_forms <- list(
        "log(serBilir)" = ~ value(log(serBilir)):CR,
        "prothrombin" = ~ value(prothrombin):CR
    )

    jointFit2 <- jm(CoxFit_CR, list(fm1, fm2), time_var = "year",
                    functional_forms = CR_forms)

    save(jointFit1, jointFit2,
         file = "C:/Users/drizo/OneDrive/Desktop/predict_simpleExample.RData")

    load("C:/Users/drizo/OneDrive/Desktop/predict_simpleExample.RData")
    load("C:/Users/drizo/OneDrive/Desktop/predict_CRExample.RData")
    data("pbc2", package = "JM")
    data("pbc2.id", package = "JM")

    pbc2.idCR <- JMbayes2::crisk_setup(pbc2.id, statusVar = "status", censLevel = "alive",
                             nameStrata = "CR")

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
    source("./R/predict_funs.R")
    source("./R/create_Wlong_mats.R")
    source("./R/help_functions.R")
    source("./R/predict_funs.R")
    Rcpp::sourceCpp('src/mcmc_fit.cpp')

    source("./Development/CI/prepare_data.R")
    source("./Development/CI/causal_effects_fun.R")
    newdataL <- pbc2[pbc2$id %in% c(81), ]
    newdataL$status2 <- 0
    # The data.frame 'newdataE_withIE' contains the event information in the
    # presence of the IE.
    pbc2_CR$serBilir <- 0.1
    newdataE_withIE <- pbc2_CR[pbc2_CR$id == 81, ]
    newdataE_withIE$event <- 0
    # The data.frame 'newdataE_withoutIE' contains the event information in the
    # absence of the IE.
    newdataE_withoutIE <- newdataE_withIE[c(1, 3), ]

    t0 <- 5

    newdataL_i <- newdataL[newdataL$year <= t0, ]
    # In the event data.frame without the IE we set the stop time at t0
    newdataE_withoutIE_i <- newdataE_withoutIE
    newdataE_withoutIE_i$stop <- t0
    # In the event data.frame with the IE we set that the IE occurs a bit after t0
    # and a bit afterward is the last time the patient was available
    newdataE_withIE_i <- newdataE_withIE
    newdataE_withIE_i$stop[c(2, 4)] <- t0

    # We calculate the predictions using the two datasets
    newdata_withIE_i <- list(newdataL = newdataL_i, newdataE = newdataE_withIE_i)
    newdata_withoutIE_i <- list(newdataL = newdataL_i, newdataE = newdataE_withoutIE_i)

}


object <- jointFit
control = NULL
#ND <- pbc2[pbc2$id %in% c(2, 3), ]
#ND$id <- factor(ND$id)
#ND <- ND[ND$year < 1, ]
#ND$status2 <- 0
#ND$years <- 1 #with(ND, ave(year, id, FUN = function (x) max(x, na.rm = T)))
#ND. <- pbc2.idCR[pbc2.idCR$id %in% c(2, 3), ]
#ND.$id <- factor(ND.$id)
#ND.$status2 <- 0
#ND.$years <- 1
#newdata = list(newdataL = ND, newdataE = ND.)
newdata = newdata_withoutIE_i
newdata2 = newdata_withIE_i
times = c(6, 7, 8)
process = "event"
type_pred = "response"
type = "subject_specific"
parallel = "snow"
level = 0.95; return_newdata = FALSE
n_samples = 500L; n_mcmc = 55L; cores = NULL; seed = 123L; use_Y = TRUE
all_times = FALSE; times_per_id = FALSE
return_mcmc = FALSE

#############################################################
#############################################################

x = predLong
x2 = predEvent
subject = 1; outcomes = 1:2
fun_long = NULL; fun_event = NULL
CI_long = TRUE; CI_event = TRUE
xlab = "Follow-up Time"; ylab_long = NULL
ylab_event = "Cumulative Risk"; main = ""
lwd_long = 2; lwd_event = 2
ylim_long_outcome_range = FALSE
col_line_long = "#0000FF"
col_line_event = c("#FF0000", "#03BF3D", "#8000FF")
pch_points = 16; col_points = "blue"; cex_points = 1
fill_CI_long = "#0000FF4D"
fill_CI_event = c("#FF00004D", "#03BF3D4D", "#8000FF4D")
cex_xlab = 1; cex_ylab_long = 1; cex_ylab_event = 1
cex_main = 1; cex_axis = 1; col_axis = "black"
pos_ylab_long = c(0.1, 2, 0.08); bg = "white"


#############################################################
#############################################################

object = jointFit1
newdata = pbc2
Tstart = 5
Thoriz = NULL
Dt = 2
type_weights = "IPCW"

tvAUC(jointFit1, newdata = pbc2, Tstart = 3, Dt = 2)
tvAUC(jointFit1, newdata = pbc2, Tstart = 3, Dt = 2, type_weights = "I")

xx1 <- tvROC(jointFit1, newdata = pbc2, Tstart = 3, Dt = 2)
xx2 <- tvROC(jointFit1, newdata = pbc2, Tstart = 3, Dt = 2, type_weights = "I")

plot(xx1)
plot(xx2)


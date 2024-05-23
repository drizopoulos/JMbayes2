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
}


object <- jointFit2
ND <- pbc2[pbc2$id %in% c(2, 3), ]
ND$id <- factor(ND$id)
ND <- ND[ND$year < 1, ]
ND$status2 <- 0
ND$years <- 1 #with(ND, ave(year, id, FUN = function (x) max(x, na.rm = T)))
ND. <- pbc2.idCR[pbc2.idCR$id %in% c(2, 3), ]
ND.$id <- factor(ND.$id)
ND.$status2 <- 0
ND.$years <- 1
newdata = list(newdataL = ND, newdataE = ND.)
newdata2 = NULL
times = c(2, 3)
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


##########################################################################################
# Author: D. Rizopoulos                                                                  #
# Aim: Test function jm()                                                                #
##########################################################################################

library("survival")
library("nlme")
library("GLMMadaptive")
library("splines")
#library("Formula")
data("pbc2", package = "JM")
data("pbc2.id", package = "JM")
source(file.path(getwd(), "R/jm.R"))
source(file.path(getwd(), "R/help_functions.R"))
source(file.path(getwd(), "Development/jm/R_to_Cpp.R"))
source(file.path(getwd(), "Development/jm/PBC_data.R"))

##########################################################################################
##########################################################################################

fm1 <- lme(log(serBilir) ~ year * sex + I(year^2) + age + prothrombin,
           data = pbc2, random = ~ year | id)
fm2 <- lme(serChol ~ ns(year, 3) + sex + age, data = pbc2, random = ~ year | id,
           na.action = na.exclude)
fm3 <- mixed_model(hepatomegaly ~ sex + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())
fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())

CoxFit <- coxph(Surv(years, status2) ~ 1,
                data = pbc2.id, model = TRUE)

##########################################################################################
##########################################################################################

Surv_object = CoxFit
Mixed_objects = list(fm1, fm2, fm3, fm4)
time_var = "year"
functional_forms = list("log(serBilir)" = ~ value(log(serBilir)) + slope(log(serBilir)) +
                            value(log(serBilir)):sex,
                        "serChol" = ~ value(serChol) + slope(serChol),
                        #"hepatomegaly" = ~ value(hepatomegaly),
                        "ascites" = ~ value(ascites) + area(ascites))
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL







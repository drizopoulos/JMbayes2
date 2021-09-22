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


#
pbc2$prothrombin[pbc2$id == levels(pbc2$id)[1L]] <- NA
pbc2$prothrombin[pbc2$id == levels(pbc2$id)[2L]] <- NA

fm1 <- lme(log(serBilir) ~ ns(year, 2, B = c(0, 15)) * sex, data = pbc2,
           random = ~ ns(year, 2, B = c(0, 15)) | id)

fm2 <- lme(prothrombin ~ ns(year, 2, B = c(0, 15)) * sex, data = pbc2,
           random = ~ ns(year, 2, B = c(0, 15)) | id,
           na.action = na.exclude, control = lmeControl(opt = "optim"))

fm3 <- mixed_model(ascites ~ year * sex, data = pbc2,
                   random = ~ 1 | id, family = binomial())
Mixed <- list(fm1, fm2, fm3)
Cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

#system.time(obj <- jm(Cox, Mixed, time_var = "year"))

#summary(obj)
#traceplot(obj, "alphas")
#gelman_diag(obj)

fForms <- ~ value(log(serBilir)):I(1 * (sex == 'male')) + value(prothrombin) + value(ascites)
fForms <- ~ value(log(serBilir)) + vsquare(value(log(serBilir))) +
    vcubic(value(log(serBilir))) + slope(log(serBilir)) +
    value(log(serBilir)):slope(log(serBilir))
fForms <- ~ value(log(serBilir)) + slope(log(serBilir))

# system.time(obj <- jm(Cox, Mixed, time_var = "year", functional_forms = fForms))


Surv_object = Cox
Mixed_objects = Mixed
time_var = 'year'
functional_forms = fForms
recurrent = FALSE
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
                   slope(log(serBilir)):sex,
               "serChol" = ~ value(serChol) + slope(serChol),
               "hepatomegaly" = ~ vexpit(value(hepatomegaly)) + sex,
               "ascites" = ~ Dexpit(value(ascites)):slope(ascites) + area(ascites))

test <- jm(CoxFit, list(fm1, fm2, fm3, fm4), time_var = "year",
           functional_forms = fForms)

tt <- jm(CoxFit, fm1, time_var = "year",
         functional_forms = ~ value(log(serBilir)) + slope(log(serBilir)))

####

Surv_object = CoxFit
Mixed_objects = list(fm1, fm2, fm3, fm4)
time_var = "year"
functional_forms = fForms
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL


################################################################################
################################################################################


fm <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
CoxFit <- coxph(Surv(years, status2) ~ 1, data = pbc2.id)
ff <- list("log(serBilir)" = ~ area(log(serBilir)))
test <- jm(CoxFit, fm, time_var = "year")
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

fm <- lme(sqrt(CD4) ~ obstime * drug, data = aids, random = ~ obstime | patient)
gm <- coxph(Surv(Time, death) ~ drug, data = aids.id)
jmFit <- jm(gm, fm, time_var = "obstime",
            functional_forms = ~ value(sqrt(CD4)) + square(value(sqrt(CD4))))
summary(jmFit)

traceplot(jmFit)

Surv_object = gm
Mixed_objects = fm
time_var = "obstime"
functional_forms = ~ value(sqrt(CD4)) + square(value(sqrt(CD4)))
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con




################################################################################

pbc2$serCholD <- as.numeric(pbc2$serChol > 250)
fm <- mixed_model(serCholD ~ year + age, data = pbc2,
                  random = ~ year | id, family = binomial())
Cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

jFit1 <- jm(Cox, fm, time_var = "year", n_iter = 6500L, n_burnin = 2000L,
            functional_forms = ~ value(serCholD))

jFit2 <- jm(Cox, fm, time_var = "year", n_iter = 6500L, n_burnin = 2000L,
            functional_forms = ~ vexpit(value(serCholD)))

jFit3 <- jm(Cox, fm, time_var = "year", n_iter = 6500L, n_burnin = 2000L,
            functional_forms = ~ vexpit(value(serCholD)) + Dexpit(slope(serCholD)))

summary(jFit1)
summary(jFit2)
summary(jFit3)

traceplot(jFit1, "alphas")
traceplot(jFit2, "alphas")
traceplot(jFit3, "alphas")


fForms <- list("ascites" = ~ vexpit(value(ascites)))
fForms <- list("ascites" = ~ vexpit(value(ascites)):slope(ascites) +
                   value(ascites)*sex + area(ascites) + slope(ascites))

pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
fm2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)
fm3 <- mixed_model(ascites ~ year + sex, data = pbc2,
                   random = ~ year | id, family = binomial())

fForms <- ~ log(value(log(serBilir)))

jointFit3 <- update(jointFit2, functional_forms = fForms)
summary(jointFit3)


Surv_object = CoxFit
Mixed_objects = list(fm1, fm2, fm3)
time_var = 'year'
functional_forms = fForms
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con
control$n_chains = 1


fm <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
jm1 <- jm(CoxFit, fm, time_var = "year", cores = 1L)
jm2 <- jm(CoxFit, fm, time_var = "year", cores = 2L)
jm3 <- jm(CoxFit, fm, time_var = "year", cores = 3L)

summary(jm1)
summary(jm2)
summary(jm3)

length(jm1$mcmc$bs_gammas)
length(jm2$mcmc$bs_gammas)
length(jm3$mcmc$bs_gammas)





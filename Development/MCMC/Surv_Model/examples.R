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

CoxFit <- coxph(Surv(years, status2) ~ age,
                data = pbc2.id, model = TRUE)

CoxFit <- survreg(Surv(years, yearsU, status3, type = "interval") ~ 1,
                  data = pbc2.id, model = TRUE)

fForms <- list("log(serBilir)" = ~ value(log(serBilir)) + slope(log(serBilir)) +
                   value(log(serBilir)):sex,
               "serChol" = ~ value(serChol) + slope(serChol),
               #"hepatomegaly" = ~ value(hepatomegaly),
               "ascites" = ~ value(ascites) + area(ascites))

test <- jm(CoxFit, list(fm1, fm2, fm3, fm4), time_var = "year",
           functional_forms = fForms)


# parameter values
betas <- test$initial_values$betas
b <- test$initial_values$b
gammas <- test$initial_values$gammas
bs_gammas <- test$initial_values$bs_gammas
alphas <- test$initial_values$alphas

# outcome vectors and design matrices
idT <- test$model_data$idT
Time_right <- test$model_data$Time_right
Time_left <- test$model_data$Time_left
Time_start <- test$model_data$Time_start
delta <- test$model_data$delta
which_event <- test$model_data$which_event
which_right <- test$model_data$which_right
which_left <- test$model_data$which_left
which_interval <- test$model_data$which_interval
W0_H <- test$model_data$W0_H
W_H <- test$model_data$W_H
X_H <- test$model_data$X_H
Z_H <- test$model_data$Z_H
U_H <- test$model_data$U_H
W0_h <- test$model_data$W0_h
W_h <- test$model_data$W_h
X_h <- test$model_data$X_h
Z_h <- test$model_data$Z_h
U_h <- test$model_data$U_h
W0_H2 <- test$model_data$W0_H2
W_H2 <- test$model_data$W_H2
X_H2 <- test$model_data$X_H2
Z_H2 <- test$model_data$Z_H2
U_H2 <- test$model_data$U_H2
log_Pwk <- test$model_data$log_Pwk
log_Pwk2 <- test$model_data$log_Pwk2





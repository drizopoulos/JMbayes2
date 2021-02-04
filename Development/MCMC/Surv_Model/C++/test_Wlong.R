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

fm1 <- lme(log(serBilir) ~ year * sex + I(year^2) + age + prothrombin,
           data = pbc2, random = ~ year | id)
fm2 <- lme(serChol ~ ns(year, 3) + sex + age, data = pbc2, random = ~ year | id,
           na.action = na.exclude)
fm3 <- mixed_model(hepatomegaly ~ sex + age + year, data = pbc2,
                   random = ~ year | id, family = binomial())
fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ year | id, family = binomial())

CoxFit <- coxph(Surv(years, status2) ~ age, data = pbc2.id)


fForms <- list("log(serBilir)" = ~ value(log(serBilir)) + slope(log(serBilir)) +
                   value(log(serBilir)):sex,
               "serChol" = ~ value(serChol) + slope(serChol),
               "hepatomegaly" = ~ vexpit(value(hepatomegaly)) + sex,
               "ascites" = ~ dexpit(value(ascites)):slope(ascites) + area(ascites))


Surv_object = CoxFit
Mixed_objects = list(fm1, fm2, fm3, fm4)
time_var = "year"
functional_forms = fForms
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con
#
betas <- initial_values$betas
b <- initial_values$b
gammas <- initial_values$gammas
bs_gammas <- initial_values$bs_gammas
alphas <- initial_values$alphas
# outcome vectors and design matrices
Time_right <- model_data$Time_right
Time_left <- model_data$Time_left
Time_start <- model_data$Time_start
delta <- model_data$delta
which_event <- model_data$which_event
which_right <- model_data$which_right
which_left <- model_data$which_left
which_interval <- model_data$which_interval
W0_H <- model_data$W0_H
W_H <- model_data$W_H
X_H <- model_data$X_H
Z_H <- model_data$Z_H
U_H <- model_data$U_H
W0_h <- model_data$W0_h
W_h <- model_data$W_h
X_h <- model_data$X_h
Z_h <- model_data$Z_h
U_h <- model_data$U_h
W0_H2 <- model_data$W0_H2
W_H2 <- model_data$W_H2
X_H2 <- model_data$X_H2
Z_H2 <- model_data$Z_H2
U_H2 <- model_data$U_H2
# other information
n <- model_data$n
idT <- model_data$idT
log_Pwk <- model_data$log_Pwk
log_Pwk2 <- model_data$log_Pwk2
FunForms_per_outcome <- model_info$FunForms_per_outcome
Funs_FunForms <- model_info$Funs_FunForms
# id_H is used to repeat the random effects of each subject GK_k times
id_H <- rep(list(rep(unclass(idT), each = control$GK_k)), length(X_H))
# this is the linear predictor for the longitudinal outcomes evaluated at the
# Gauss-Kronrod quadrature points
eta_H <- linpred_surv(X_H, betas, Z_H, b, id_H)
# Wlong is the design matrix of all longitudinal outcomes according to the specified
# functional forms per outcome already multiplied with the interaction terms matrix U
Wlong_H <- create_Wlong(eta_H, FunForms_per_outcome, U_H, Funs_FunForms)


DD <- list(eta = eta_H, U = U_H, FunForms_cpp = model_info$FunForms_cpp,
           FunForms_ind = model_info$FunForms_ind,
           Funs_FunForms = model_info$Funs_FunForms)

test <- create_Wlong_cpp(DD)


all.equal(Wlong_H[[1]], test[[1]], check.attributes = FALSE)
all.equal(Wlong_H[[2]], test[[2]], check.attributes = FALSE)
all.equal(Wlong_H[[3]], test[[3]], check.attributes = FALSE)
all.equal(Wlong_H[[4]], test[[4]], check.attributes = FALSE)

library("rbenchmark")
benchmark(R = create_Wlong(eta_H, FunForms_per_outcome, U_H, Funs_FunsForms),
          Cpp = create_Wlong_cpp(DD), replications = 10000)




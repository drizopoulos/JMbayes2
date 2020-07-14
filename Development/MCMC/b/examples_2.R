library("survival")
library("nlme")
library("GLMMadaptive")
library("splines")
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

CoxFit <- coxph(Surv(years, status2) ~ age * sex,
                data = pbc2.id, model = TRUE)

#CoxFit <- survreg(Surv(years, yearsU, status3, type = "interval") ~ 1,
#                  data = pbc2.id, model = TRUE)

fForms <- list("log(serBilir)" = ~ value(log(serBilir)) + slope(log(serBilir)) +
                 value(log(serBilir)):sex,
               "serChol" = ~ value(serChol) + slope(serChol),
               "hepatomegaly" = ~ value(hepatomegaly),
               "ascites" = ~ value(ascites) + area(ascites))

JM1 <- jm(CoxFit, list(fm1, fm2, fm3, fm4), time_var = "year",
          functional_forms = fForms)

#JM1 <- jm(CoxFit, list(fm1, fm3, fm4), time_var = "year",
#          functional_forms = fForms)

###########################################################

test <- JM1

# parameter values
betas <- test$initial_values$betas
b <- test$initial_values$b
gammas <- test$initial_values$gammas
bs_gammas <- test$initial_values$bs_gammas
alphas <- test$initial_values$alphas

# outcome vectors and design matrices
n <- test$model_data$n
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
Xhc <- test$model_data$Xhc
columns_HC <- test$model_data$columns_HC
unq_idL <- test$model_data$unq_idL
X <- test$model_data$X
Z <- test$model_data$Z
idL_lp <- test$model_data$idL_lp
y <- test$model_data$y
nY <- length(unique(test$model_data$unq_idL[[1]]))
log_sigmas <- test$initial_values$log_sigmas
idL <- test$model_data$idL
Funs <- lapply(test$model_info$families, log_dens_Funs)
mu_funs <- lapply(test$model_info$families, "[[", 'linkinv')

control <- test$control
functional_forms_per_outcome <- test$model_info$fun_forms$functional_forms_per_outcome

# id_H is used to repeat the random effects of each subject GK_k times
id_H <- lapply(X_H, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
# this is the linear predictor for the longitudinal outcomes evaluated at the
# Gauss-Kronrod quadrature points
eta_H <- linpred_surv(X_H, betas, Z_H, b, id_H)
# Wlong is the design matrix of all longitudinal outcomes according to the specified
# functional forms per outcome already multiplied with the interaction terms matrix U
Wlong_H <- create_Wlong(eta_H, functional_forms_per_outcome, U_H)
if (length(which_event)) {
  id_h <- lapply(X_h, function (x) seq_len(nrow(x[[1]])))
  eta_h <- linpred_surv(X_h, betas, Z_h, b, id_h)
  Wlong_h <- create_Wlong(eta_h, functional_forms_per_outcome, U_h)
}
if (length(which_interval)) {
  id_H2 <- lapply(X_H2, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
  eta_H2 <- linpred_surv(X_H2, betas, Z_H, b, id_H2)
  Wlong_H2 <- create_Wlong(eta_H2, functional_forms_per_outcome, U_H2)
} else {
  Wlong_H2 <- rep(list(matrix(0.0, length(Time_right), 1)), length(W_H))
}

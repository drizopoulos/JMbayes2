# Test Mstate with current version
library(JMbayes2)
library(JMbayes)
library(mstate)
library(nlme)
library(mstate) # Please use the version 0.2.7
library(devtools)
library(JM)



# LOIC IMPLEMENTATION
# Construct the 3*3 matrix of transitions:
tmat <- matrix(NA, 3, 3)
tmat[1, 2:3] <- 1:2
tmat[2, 3] <- 3
dimnames(tmat) <- list(from = c("State_0", "State_1", "State_2"),
                       to = c("State_0", "State_1", "State_2"))
tmat

covs <- "X"

data_mstate <- msprep(time = c(NA, "t_State_1", "t_State_2"),
                      status = c(NA, "State_1", "State_2"),
                      data = data_surv,
                      trans = tmat,
                      keep = covs,
                      id = "id")

data_mstate <- msprep(time = c(NA, "t_State_1", "t_State_2"),
                      status = c(NA, "State_1", "State_2"),
                      data = data_surv,
                      trans = tmat,
                      id = "id")
data_mstate

data_mstate <- expand.covs(data_mstate, covs,
                           append = TRUE, longnames = FALSE)


lmeFit <- lme(fixed = Y ~ (times + I((1 + times)^(-1.2) - 1)) * X,
              data = data_long,
              random = ~ (times + I((1 + times)^(-1.2) - 1)) | id,
              method = "REML",
              control = list(opt = "optim"))

coxFit <- coxph(Surv(Tstart, Tstop, status) ~ X.1 + X.2 + X.3 + strata(trans),
                data = data_mstate, method = "breslow", x = TRUE, model = TRUE)





dForm <- list(fixed = ~ 1 + I((-1.2) * ((1 + times)^(-2.2))) +
                  X + I((-1.2) * ((1 + times)^(-2.2))):X,
              indFixed = c(2:3 ,5:6),
              random = ~ 1 + I((-1.2) * ((1 + times)^(-2.2))),
              indRandom = 2:3)


jointFit_1step_GHk3 <-
    JMstateModel(lmeObject = lmeFit,
                 survObject = coxFit,
                 timeVar = "times",
                 parameterization = "both",
                 method = "spline-PH-aGH",
                 interFact = list(value = ~ strata(trans) - 1,
                                  slope = ~ strata(trans) - 1,
                                  data = data_mstate),
                 derivForm = dForm,
                 Mstate = TRUE,
                 data.Mstate = data_mstate,
                 ID.Mstate = "id",
                 control = list(GHk = 3, lng.in.kn = 1),
                 verbose = TRUE)

# JMbayes implementation
mixed_model <- mvglmer(list(Y ~ (times + I((1 + times)^(-1.2) - 1)) * X + ((times + I((1 + times)^(-1.2) - 1)) | id)),
                       families = list(gaussian),
                       data = data_long)

coxFit2 <- coxph(Surv(Tstart, Tstop, status) ~ X.1 + X.2 + X.3 + strata(trans) + cluster(id),
                 data = data_mstate, method = "breslow", x = TRUE, model = TRUE)

interacts <- list("Y" = ~ strata(trans) - 1)

jm_mstate_model <- mvJointModelBayes(mixed_model, coxFit2, timeVar = 'times',
                                     Interactions = interacts,
                                     multiState = TRUE,
                                     data_MultiState = data_mstate,
                                     idVar_MultiState = 'id')

# JMbayes2 test implementation
ms_forms <- list(
    "Y" = ~ value(Y):trans
)

jFit_ms <- jm(coxFit2, list(lmeFit), time_var = "times",
              functional_forms = ms_forms)

summary(jFit_CR)

###############################################################


lmeFit_DR <- lme(fixed = Y ~ times * X, data = data_long,
                 random = ~ times | id, control = list(opt = "optim"))

data_mstate$trans <- factor(data_mstate$trans)
CoxFit_DR <- coxph(Surv(Tstart, Tstop, status) ~ X * strata(trans),
                   data = data_mstate)


jmFit_DR <- jm(CoxFit_DR, lmeFit_DR, time_var = "times",
               functional_forms = ~ value(Y):trans)

summary(jmFit_DR)
traceplot(jmFit_DR)


################# COMPARE wih JMstateModel
CoxFit_DR_adj <- coxph(Surv(Tstart, Tstop, status) ~ X * strata(trans),
                       data = data_mstate, method = "breslow", x = TRUE, model = TRUE)


jointFit_1step_GHk3_adj <-
  JMstateModel(lmeObject = lmeFit_DR,
               survObject = CoxFit_DR_adj,
               timeVar = "times",
               parameterization = "value",
               method = "spline-PH-aGH",
               interFact = list(value = ~ strata(trans) - 1,
                                data = data_mstate),
               Mstate = TRUE,
               data.Mstate = data_mstate,
               ID.Mstate = "id",
               control = list(GHk = 3, lng.in.kn = 1),
               verbose = TRUE)

summary(jointFit_1step_GHk3_adj)

mixed_model_adj <- mvglmer(list(Y ~ times * X + (times | id)),
                       families = list(gaussian),
                       data = data_long)

CoxFit_DR_adj <- coxph(Surv(Tstart, Tstop, status) ~ X * strata(trans) + cluster(id),
                   data = data_mstate, model = TRUE, x = TRUE)

interacts_adj <- list("Y" = ~ strata(trans) - 1)

jm_mstate_model_adj <- mvJointModelBayes(mixed_model_adj, CoxFit_DR_adj, timeVar = 'times',
                                     Interactions = interacts_adj,
                                     multiState = TRUE,
                                     data_MultiState = data_mstate,
                                     idVar_MultiState = 'id')


jointFit_1step_GHk3_adj$coefficients$alpha
jmFit_DR$statistics$Mean$alphas
jm_mstate_model_adj$statistics$postMeans$alphas

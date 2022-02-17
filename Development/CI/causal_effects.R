library("JMbayes2")
library("lattice")
source("./Development/CI/prepare_data.R")
source("./Development/CI/causal_effects_fun.R")


####################################
# Causal Effects from Joint Models #
####################################

# We fit the joint model of interest.

# In the longitudinal sub-model we postulate that the trajectory changes
# after the occurrence of the intermediate event IE. Note that the formulation
# of the model allows for a 'sudden' drop/jump of the subject-specific profiles,
# i.e., the intercept and slope change.
lmeFit <- lme(log(serBilir) ~ year + IE + IE:year, data = pbc2,
              random = ~ year + IE + IE:year | id,
              control = lmeControl(opt = "optim"))

# In the event process sub-model we fit cause-specific hazards functions, also
# including IE as a time-varying covariate
CoxFit <- coxph(Surv(start, stop, event) ~ strata(CR) * (age + sex + IE),
                data = pbc2_CR)

# In the functional forms we specify that the effect of the longitudinal
# biomarker is different before and after the IE status for death, but the
# biomarker has no effect on the risk of transplantation
dummy <- function (f, lvl) as.numeric(f == lvl)
fForms <- ~ value(log(serBilir)):dummy(CR, "transplanted") +
    value(log(serBilir)):dummy(CR, "dead"):dummy(IE, 0) +
    value(log(serBilir)):dummy(CR, "dead"):dummy(IE, 1)

# To specify that the biomarker has no effect on the risk of transplantation, we
# constraint the corresponding association parameter alpha via the prior (i.e.,
# we assume prior mean zero and prior precision 20000)
jointFit <- jm(CoxFit, lmeFit, time_var = "year", functional_forms = fForms,
               priors = list(Tau_alphas = list(diag(20000, 1), diag(1), diag(1))),
               n_iter = 15000L, n_burnin = 5000L)

summary(jointFit)

#-----------------------------
# Conditional Causal Effects -
#-----------------------------

# We take as an example Patient 81 from the PBC dataset.
# We calculate the cumulative incidence probabilities in the presence and
# absence of the intermediate event IE.
# The data.frame 'newdataL' contains the longitudinal measurements
newdataL <- pbc2[pbc2$id %in% c(81), ]
newdataL$status2 <- 0
# The data.frame 'newdataE_withIE' contains the event information in the
# presence of the IE.
newdataE_withIE <- pbc2_CR[pbc2_CR$id %in% c(81), ]
newdataE_withIE$event <- 0
# The data.frame 'newdataE_withoutIE' contains the event information in the
# absence of the IE.
newdataE_withoutIE <- newdataE_withIE[c(1, 3), ]

# We specify the follow-up times t0 at which we want to calculate the
# conditional causal effect, and the length of the medically-relevant time
# window Delta_t
t0 <- c(3, 5, 7 , 9)
Delta_t <- 2

conditional_causal_effect <-
    matrix(0.0, length(t0), 3,
           dimnames = list(paste("t0 =", t0),
                           c("Effect", "95%_low", "95%_upp")))

for (i in seq_along(t0)) {
    # Some data management steps first:
    # In the longitudinal data.frame we keep the measurements before t0
    newdataL_i <- newdataL[newdataL$year <= t0[i], ]
    # In the event data.frame without the IE we set the stop time at t0
    newdataE_withoutIE_i <- newdataE_withoutIE
    newdataE_withoutIE_i$stop <- t0[i]
    # In the event data.frame with the IE we set that the IE occurs a bit after t0
    # and a bit afterward is the last time the patient was available
    newdataE_withIE_i <- newdataE_withoutIE_i
    newdataE_withIE_i$IE <- 1

    # We calculate the predictions using the two datasets
    newdata_withIE_i <- list(newdataL = newdataL_i, newdataE = newdataE_withIE_i)
    newdata_withoutIE_i <- list(newdataL = newdataL_i, newdataE = newdataE_withoutIE_i)

    predSurv_withoutIE <- predict(jointFit, newdata = newdata_withoutIE_i,
                                  process = "event", times = t0[i] + Delta_t,
                                  return_mcmc = TRUE)
    predSurv_withIE <- predict(jointFit, newdata = newdata_withoutIE_i,
                               newdata2 = newdata_withIE_i, return_mcmc = TRUE,
                               process = "event", times = t0[i] + Delta_t)

    # The conditional causal effect
    conditional_causal_effect[i, 1] <-
        predSurv_withIE$pred[2] - predSurv_withoutIE$pred[2]


    # test <- causal_effects(jointFit, newdataL_i, newdataE_withoutIE_i,
    #                       newdataE_withIE_i, t0 = t0[i], Dt = Delta_t)

    # 95% CI
    samples <- predSurv_withIE$mcmc[2, ] - predSurv_withoutIE$mcmc[2, ]
    conditional_causal_effect[i, 2:3] <-
        quantile(samples, probs = c(0.025, 0.975))
}

conditional_causal_effect

#--------------------------
# Marginal Causal Effects -
#--------------------------

# To calculate the marginal causal effect, we will average over the subjects
# at risk,

# We specify the follow-up time t0 at which we want to calculate the
# marginal causal effect, and the length of the medically-relevant time
# window Delta_t

# we get the data of the subjects at risk at t0 and with two datasets for the
# event outcome, one in which the IE is set at zero and on in which it is set
# at one
Data <- get_data(pbc2, pbc2_CR, t0 = 3, Dt = 2, object = jointFit, IE_var = "IE")


# we calculate the effects
effects <- causal_effects(jointFit, Data$newdataL, Data$newdataE, Data$newdataE2,
                          t0 = 3, Dt = 2, extra_objects = "dummy", B = 5,
                          calculate_CI = TRUE)






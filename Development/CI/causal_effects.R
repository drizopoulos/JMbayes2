library("JMbayes2")
library("lattice")
source("./Development/CI/prepare_data.R")
source("./Development/CI/causal_effects_fun.R")

####################################
# Causal Effects from Joint Models #
####################################

# We fit the joint model of interest.

# In the longitudinal sub-model, we postulate that the trajectory changes
# after the occurrence of the intermediate event IE. Note that the formulation
# of the model allows for a 'sudden' drop/jump of the subject-specific profiles,
# i.e., the intercept and slope change.
lmeFit <- lme(log(serBilir) ~ year + IE + IE:year, data = pbc2,
              random = ~ year + IE + IE:year | id,
              control = lmeControl(opt = "optim"))

# In the event process sub-model, we fit cause-specific hazards functions, also
# including IE as a time-varying covariate
CoxFit <- coxph(Surv(start, stop, event) ~ strata(CR) * (age + sex + IE),
                data = pbc2_CR)

# In the functional forms, we specify that the effect of the longitudinal
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
# The data.frame 'newdataL' contains the longitudinal measurements of this patient
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
    newdataE_withIE_i <- newdataE_withIE
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
    #                      newdataE_withIE_i, t0 = t0[i], Dt = Delta_t)

    # 95% CI
    # the samples object contains the MCMC samples for the conditional
    # causal effect
    samples <- predSurv_withIE$mcmc[2, ] - predSurv_withoutIE$mcmc[2, ]

    M <- 20L
    conditional_causal_effect. <- numeric(M)
    for (m in seq_len(M)) {
        set.seed(2022L + m)
        # we want to calculate the extra variance due to the fact that we
        # condition on Y_i(t). We first simulate new Y_i(t) data and we put
        # them in the corresponding datasets.
        mu <- predict(jointFit, newdata = newdata_withoutIE_i)$preds[[1]]
        sigma <- jointFit$statistics$Mean$sigmas
        preds <- exp(rnorm(length(mu), mu, sigma))
        newdata_withoutIE_i. <- newdata_withoutIE_i
        newdata_withoutIE_i.$newdataL$serBilir <- preds
        newdata_withIE_i. <- newdata_withIE_i
        newdata_withIE_i.$newdataL$serBilir <- preds
        # we calculate the conditional causal effect using the new simulated data
        predSurv_withoutIE. <- predict(jointFit, newdata = newdata_withoutIE_i.,
                                       process = "event", times = t0[i] + Delta_t)
        predSurv_withIE. <- predict(jointFit, newdata = newdata_withoutIE_i.,
                                    newdata2 = newdata_withIE_i.,
                                    process = "event", times = t0[i] + Delta_t)

        # The conditional causal effect
        conditional_causal_effect.[m] <-
            predSurv_withIE.$pred[2] - predSurv_withoutIE.$pred[2]
    }
    total_variance <- var(conditional_causal_effect.) + var(samples)

    conditional_causal_effect[i, 2:3] <-
        c(conditional_causal_effect[i, 1] - 1.96 * sqrt(total_variance),
          conditional_causal_effect[i, 1] + 1.96 * sqrt(total_variance))
}

conditional_causal_effect

#--------------------------
# Marginal Causal Effects -
#--------------------------

# To calculate the marginal causal effect, we will average over the
# respective group subjects

# We specify the follow-up time t0 at which we want to calculate the
# marginal causal effect, and the length of the medically-relevant time
# window Dt

# We get the data of the subjects at risk at t0 and who did not have the IE yet;
# we also create two datasets for the event outcome, one in which the IE is set
# at zero and on in which it is set at one. The second one is used in the
# 'newdata2' argument of the predict() method.
Data <- get_data(pbc2, pbc2_CR, t0 = 3, Dt = 2, object = jointFit,
                 IE_var = "IE", IE_time = "S")

# we calculate the effects
marginal_causal_effect <-
    causal_effects(jointFit, Data$newdataL, Data$newdataL2,
                   Data$newdataE, Data$newdataE2,
                   t0 = 3, Dt = 2, extra_objects = "dummy", B = 5,
                   calculate_CI = TRUE)

marginal_causal_effect

#--------------------------------------
# Marginal-Conditional Causal Effects -
#--------------------------------------

# To calculate the marginal-conditional causal effect, we will average over
# the respective group subjects

# We specify the follow-up time t0 at which we want to calculate the
# marginal causal effect, and the length of the medically-relevant time
# window Dt

# We get the data of the subjects at risk at t0 and who did not have the IE yet;
# we also create two datasets for the event outcome, one in which the IE is set
# at zero and on in which it is set at one. The second one is used in the
# 'newdata2' argument of the predict() method.
Data <- get_data(pbc2, pbc2_CR, t0 = 3, Dt = 2, object = jointFit,
                 IE_var = "IE", IE_time = "S")

# This is the same as above for the marginal effect, but we now want to put
# the extra restriction that we want to keep the subjects who had their last
# serum bilirubin measurement above the value of 2.

# First, we find the last serum bilirubin measurement per subject
last <- with(Data$newdataL, tapply(serBilir, id, tail, n = 1L))
# we find id of the subjects with this measurement above 3
ids <- names(last)[last > 2]
# we keep only these subjects in the respective datasets in Data
Data$newdataL <- Data$newdataL[Data$newdataL$id %in% ids, ]
Data$newdataL2 <- Data$newdataL2[Data$newdataL2$id %in% ids, ]
Data$newdataE <- Data$newdataE[Data$newdataE$id %in% ids, ]
Data$newdataE2 <- Data$newdataE2[Data$newdataE2$id %in% ids, ]

# we calculate the effects
mc_causal_effect <-
    causal_effects(jointFit, Data$newdataL, Data$newdataL2,
                   Data$newdataE, Data$newdataE2,
                   t0 = 3, Dt = 2, extra_objects = "dummy", B = 5,
                   calculate_CI = TRUE)

mc_causal_effect

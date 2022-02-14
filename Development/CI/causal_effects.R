library("JMbayes2")
library("lattice")
source("./Development/CI/prepare_data.R")

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
t0 <- 3
Delta_t <- 2

# We are going the create a function such that we can later calculate the variance

get_marginal_effect <- function (Data_Long, Data_Event, t0, Delta_t, object) {
    # We need the subset of patients who were event-free at t0 and who did not have
    # the intermediate event yet;
    time_var <- object$model_info$var_names$time_var
    # first, we keep the longitudinal measurements before t0 for patients who
    # were at risk at t0
    R_long <- Data_Long[Data_Long[[time_var]] <= t0 & Data_Long$years > t0, ]
    R_long$id <- factor(R_long$id)
    # then we want to keep only the patients who did not have the IE yet
    keep <- function (x) rep(all(x == 0), length(x))
    no_IEs_before_t0 <- with(R_long, ave(IE, id, FUN = keep))
    R_long <- R_long[as.logical(no_IEs_before_t0), ]

    # we keep the same patients in the event times database
    R_event <- Data_Event[Data_Event$id %in% unique(R_long$id), ]
    R_event$id <- factor(R_event$id)
    # there are some patients who had an IE after t0; we need
    # to remove these lines from R_event
    R_event <- R_event[R_event$IE == 0, ]
    # we set the stop time to t0
    R_event$stop <- t0
    # we set the event to zero
    R_event$event <- 0

    ###########

    # we calculate the CIFs with IE
    newdata_withoutIE <- list(newdataL = R_long, newdataE = R_event)

    CIF_withoutIE <- predict(object, newdata = newdata_withoutIE,
                             process = "event", times = t0 + Delta_t,
                             return_mcmc = TRUE)

    # we create a copy of R_event and we set IE to one, and we calculate
    # the CIF with an IE after t0 (using 'newdata2')
    R_event2 <- R_event
    R_event2$IE <- 1
    newdata_withIE <- list(newdataL = R_long, newdataE = R_event2)

    CIF_withIE <- predict(object, newdata = newdata_withoutIE,
                          newdata2 = newdata_withIE,
                          process = "event", times = t0 + Delta_t,
                          return_mcmc = TRUE)

    # the marginal effect is the mean over the conditional effects
    mean(CIF_withIE$pred[CIF_withIE$times > t0] -
             CIF_withoutIE$pred[CIF_withoutIE$times > t0])
}

make_bootSample <- function (Data_Long, Data_Event, id_var = "id", seed = 1L) {
    if (!exists(".Random.seed", envir = .GlobalEnv))
        runif(1)
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
    ids <- Data_Long[[id_var]]
    unq_ids <- unique(ids)
    ids <- factor(ids, levels = unq_ids)
    set.seed(seed)
    new_ids <- sample(unq_ids, replace = TRUE)
    new_Data_Long <- new_Data_Event <- vector("list", length(unq_ids))
    for (i in seq_along(unq_ids)) {
        keep <- Data_Long[[id_var]] == new_ids[i]
        dataL_i <- Data_Long[keep, ]
        dataL_i[[id_var]] <- i
        new_Data_Long[[i]] <- dataL_i
        ##
        keep <- Data_Event[[id_var]] == new_ids[i]
        dataE_i <- Data_Event[keep, ]
        dataE_i[[id_var]] <- i
        new_Data_Event[[i]] <- dataE_i
    }
    list(Data_Long = do.call("rbind", new_Data_Long),
         Data_Event = do.call("rbind", new_Data_Event))
}


# the marginal effect in the original data
marginal_effect <- get_marginal_effect(pbc2, pbc2_CR, t0, Delta_t, jointFit)


# the marginal effect in M Bootstrap samples
M <- 10
Marginal_Effects <- numeric(M)
for (m in seq_len(M)) {
    Data_m <- make_bootSample(pbc2, pbc2_CR, seed = m)
    meffect <- get_marginal_effect(Data_m$Data_Long, Data_m$Data_Event,
                                   t0, Delta_t, jointFit)
    Marginal_Effects[m] <- meffect
}

var(Marginal_Effects)





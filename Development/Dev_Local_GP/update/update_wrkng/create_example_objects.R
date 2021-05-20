library('JMbayes2')
?jm

# [1] Fit the mixed model using lme().
fm1 <- lme(fixed = log(serBilir) ~ year * sex + I(year^2) +
             age + prothrombin, random =  ~ year | id, data = pbc2)

# [2] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [3] The basic joint model is fitted using a call to jm() i.e.,
joint_model_fit_1 <- jm(fCox1, fm1, time_var = "year", control = list('save_random_effects' = FALSE))
joint_model_fit_2 <- jm(fCox1, fm1, time_var = "year", control = list('save_random_effects' = TRUE))
joint_model_fit_3 <- jm(fCox1, fm1, time_var = "year", control = list('save_random_effects' = FALSE, 'n_chains' = 1))
joint_model_fit_4 <- jm(fCox1, fm1, time_var = "year", control = list('save_random_effects' = TRUE, 'n_chains' = 1))


save(fCox1, fm1, pbc2, pbc2.id, joint_model_fit_1, joint_model_fit_2, joint_model_fit_3, joint_model_fit_4, initial_values,
     file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Development\\Dev_Local_GP\\update\\update_wrkng\\example_objects_inits.RData')

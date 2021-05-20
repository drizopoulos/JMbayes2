remove.packages('JMbayes2')
remotes::install_github('drizopoulos/JMbayes2')

library(JMbayes2)
set.seed(543)
# [1] Fit the mixed model using lme().
fm1 <- lme(fixed = log(serBilir) ~ year * sex + I(year^2) +
             age + prothrombin, random =  ~ year | id, data = pbc2)

# [2] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [3] The basic joint model is fitted using a call to jm() i.e.,
joint_model_fit_1 <- jm(fCox1, fm1, time_var = "year")

summary(joint_model_fit_1)
traceplot(joint_model_fit_1)

save(joint_model_fit_1, file = paste0(getwd(), '/Development/Dev_Local_GP/update/update_wrkng/oldfit.RData'))
save(joint_model_fit_1_new, file = paste0(getwd(), '/Development/Dev_Local_GP/update/update_wrkng/newfit.RData'))

# two longitudinal outcomes example
# [1] Fit the mixed-effects models using lme() for continuous
# outcomes and mixed_model() for categorical outcomes.
fm1 <- lme(fixed = log(serBilir) ~ year * sex,
           random = ~ year | id, data = pbc2)

fm2 <- mixed_model(hepatomegaly ~ sex + age + year, data = pbc2,
                   random = ~ year | id, family = binomial())

# [2] Save all the fitted mixed-effects models in a list.
Mixed <- list(fm1, fm2)

# [3] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [4] The joint model is fitted using a call to jm() i.e.,
joint_model_fit_2 <- jm(fCox1, Mixed, time_var = "year")

save(fm1, fm2, Mixed, fCox1, joint_model_fit_2, 
     file =  paste0(getwd(), '/Development/Dev_Local_GP/update/update_new/mvfit.RData'))

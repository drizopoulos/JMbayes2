remove.packages('JMbayes2')

remotes::install_github('drizopoulos/JMbayes2')
library('JMbayes2')

?jm

# [1] Fit the mixed model using lme().
fm1 <- lme(fixed = log(serBilir) ~ year * sex + I(year^2) +
             age + prothrombin, random =  ~ year | id, data = pbc2)

# [2] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [3] The basic joint model is fitted using a call to jm() i.e.,
set.seed(12)
joint_model_fit_1_old <- jm(fCox1, fm1, time_var = "year")
joint_model_fit_1_new <- jm(fCox1, fm1, time_var = "year")


save(joint_model_fit_1_old, file = paste0(getwd(), '/Development/Dev_Local_GP/update/update_new/joint_model_fit_1_old.RData'))
save(joint_model_fit_1_new, file = paste0(getwd(), '/Development/Dev_Local_GP/update/update_new/joint_model_fit_1_new.RData'))

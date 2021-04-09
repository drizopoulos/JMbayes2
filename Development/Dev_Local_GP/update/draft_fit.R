library('JMbayes2')
?jm

# [1] Fit the mixed model using lme().
fm1 <- lme(fixed = log(serBilir) ~ year * sex + I(year^2) +
             age + prothrombin, random =  ~ year | id, data = pbc2)

# [2] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [3] The basic joint model is fitted using a call to jm() i.e.,
joint_model_fit_1 <- jm(fCox1, fm1, time_var = "year")
joint_model_fit_2 <- jm(fCox1, fm1, time_var = "year", control = list('n_iter' = 5000))



joint_model_new <- jm(fCox1, fm1, time_var = "year")
set.seed(12)
joint_model_new_2a <- jm(fCox1, fm1, time_var = "year", control = list(save_random_effects = FALSE))
set.seed(12)
joint_model_new_2b <- jm(fCox1, fm1, time_var = "year", control = list(save_random_effects = TRUE))

dim(joint_model_new_2a$mcmc$b$`1`)
dim(joint_model_new_2b$mcmc$b$`1`)

all.equal(joint_model_new_2a$mcmc$b$`1`[, ,1], joint_model_new_2b$mcmc$b$`1`[, , 6000])

joint_model_new_2a$mcmc$b$`1`[312, , 1] 
joint_model_new_2b$mcmc$b$`1`[312, , 6000]

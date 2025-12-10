fm1 <- lme(fixed = log(serBilir) ~ year * sex, 
           random =~ year | id, 
           data = pbc2)
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

fform1 <- list("log(serBilir)" = ~ value(log(serBilir)))
jm1 <- jm(fCox1, fm1, time_var = "year", data = pbc2.id,
          functional_forms = fform1)

fform2 <- list("log(serBilir)" = ~ change(log(serBilir), time_window = NULL, standardise = TRUE))
jm2 <- jm(fCox1, fm1, time_var = "year", data = pbc2.id, 
          functional_forms = fform2)

fform3 <- list("log(serBilir)" = ~ change(log(serBilir), time_window = NULL, standardise = FALSE))
jm3 <- jm(fCox1, fm1, time_var = "year", data = pbc2.id, 
          functional_forms = fform3)

fform4 <- list("log(serBilir)" = ~ change(log(serBilir), time_window = 1, standardise = TRUE))
jm4 <- jm(fCox1, fm1, time_var = "year", data = pbc2.id, 
          functional_forms = fform4)

fform5 <- list("log(serBilir)" = ~ change(log(serBilir), time_window = 2, standardise = FALSE))
jm5 <- jm(fCox1, fm1, time_var = "year", data = pbc2.id, 
          functional_forms = fform5)

fform6 <- list("log(serBilir)" = ~ change(log(serBilir), time_window = 2, standardise = FALSE) + area(log(serBilir), time_window = 2))
jm6 <- jm(fCox1, fm1, time_var = "year", data = pbc2.id, 
          functional_forms = fform6)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
rmarkdown::render("vignettes/Baseline_Hazard.Rmd")
rmarkdown::render("vignettes/Causal_Effects.Rmd")
rmarkdown::render("vignettes/Competing_Risks.Rmd")
rmarkdown::render("vignettes/Dynamic_Predictions.Rmd")
rmarkdown::render("vignettes/JMbayes2.Rmd")
rmarkdown::render("vignettes/Multi_State_Processes.Rmd")
rmarkdown::render("vignettes/Non_Gaussian_Mixed_Models.Rmd")
rmarkdown::render("vignettes/Recurring_Events.Rmd")
rmarkdown::render("vignettes/Time_Varying_Effects.Rmd")
rmarkdown::render("vignettes/Transformation_Functions.Rmd")
rmarkdown::render("vignettes/Super_Learning.Rmd")




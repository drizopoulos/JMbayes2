library(JMbayes2)
pbc2.id[pbc2.id$id %in% c(1, 2, 5), c("id", "years", "status")]
pbc2.idCR <- crisk_setup(pbc2.id, statusVar = "status", censLevel = "alive", 
                         nameStrata = "CR")
pbc2.idCR[pbc2.idCR$id %in% c(1, 2, 5), 
          c("id", "years", "status", "status2", "CR")]
CoxFit_CR <- coxph(Surv(years, status2) ~ (age + drug) * strata(CR),
                   data = pbc2.idCR)
fm1 <- lme(log(serBilir) ~ poly(year, 2) * drug, data = pbc2, 
           random = ~ poly(year, 2) | id)
fm2 <- lme(prothrombin ~ year * drug, data = pbc2, random = ~ year | id)


pbc2.idCR$IE1 <- pbc2.idCR$years
pbc2.idCR$IE1[pbc2.idCR$CR == "death" & pbc2.idCR$status2 == 0] <- Inf 
pbc2.idCR$IE2[pbc2.idCR$CR == "transplanted" & pbc2.idCR$status2 == 0] <- Inf 

CR_forms <- list(
  "log(serBilir)" = ~ value(log(serBilir), IE_time = "IE1"):CR,
  "prothrombin" = ~ value(prothrombin):CR
)
jFit_CR1 <- jm(CoxFit_CR, list(fm1, fm2), time_var = "year",
              functional_forms = CR_forms,
              n_iter = 25000L, n_burnin = 5000L, n_thin = 5L)
summary(jFit_CR1)

jFit_CR2 <- jm(CoxFit_CR, list(fm1, fm2), time_var = "year",
               functional_forms = NULL,
               n_iter = 25000L, n_burnin = 5000L, n_thin = 5L)
summary(jFit_CR2)

# Surv_object <- CoxFit_CR
# Mixed_objects <- list(fm1, fm2)
# time_var <- "year"
# recurrent <- FALSE
# functional_forms <- CR_forms
# 
# which_independent <- NULL
# base_hazard <- NULL
# data_Surv <- NULL 
# id_var <- NULL
# priors <- NULL
# control <- NULL
# 
# source("R/help_functions.R")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

library(JMbayes2)
fm1  <- lme(fixed = log(serBilir) ~ year + sex + age, 
            random =  ~ year | id, data = pbc2)
fm2 <- lme(prothrombin ~ year + sex, data = pbc2, random = ~ year | id)

pbc2.id$IE_time <- rep(Inf, 312)
pbc2.id$IE_time2 <- pbc2.id$years

cox <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)
ff1a <- list("log(serBilir)" = ~ value(log(serBilir)))
ff1b <- list("log(serBilir)" = ~ value(log(serBilir), IE_time = "IE_time"))
jm1a <- jm(cox, fm1, time_var = "year", functional_forms = ff1a)
jm1b <- jm(cox, fm1, time_var = "year", functional_forms = ff1b)
all(fixef(jm1a)[[1]] == fixef(jm1b)[[1]])


cox <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)
ff1a <- list("log(serBilir)" = ~ value(log(serBilir)))
ff1b <- list("log(serBilir)" = ~ value(log(serBilir), IE_time = "IE_time2"))
jm1a <- jm(cox, fm1, time_var = "year", functional_forms = ff1a)
jm1b <- jm(cox, fm1, time_var = "year", functional_forms = ff1b)
all(fixef(jm1a)[[1]] == fixef(jm1b)[[1]])


ff2a <- list("log(serBilir)" = ~ slope(log(serBilir)))
ff2b <- list("log(serBilir)" = ~ slope(log(serBilir), IE_time = "IE_time"))
jm2a <- jm(cox, fm1, time_var = "year", functional_forms = ff2a)
jm2b <- jm(cox, fm1, time_var = "year", functional_forms = ff2b)
all(fixef(jm2a)[[1]] == fixef(jm2b)[[1]])


ff3a <- list("log(serBilir)" = ~ area(log(serBilir)))
ff3b <- list("log(serBilir)" = ~ area(log(serBilir), IE_time = "IE_time"))
jm3a <- jm(cox, fm1, time_var = "year", functional_forms = ff3a)
jm3b <- jm(cox, fm1, time_var = "year", functional_forms = ff3b)
all(fixef(jm3a)[[1]] == fixef(jm3b)[[1]])


Surv_object <- cox
Mixed_objects <- fm1
time_var <- "year"
recurrent <- FALSE
functional_forms <- ff3b

which_independent <- NULL
base_hazard <- NULL
data_Surv <- NULL
id_var <- NULL
priors <- NULL
control <- NULL

source("R/help_functions.R")





ff4a <- list("log(serBilir)" = ~ value(log(serBilir)) + slope(log(serBilir)) + area(log(serBilir)))
ff4b <- list("log(serBilir)" = ~ value(log(serBilir), IE_time = "IE_time") + slope(log(serBilir), IE_time = "IE_time") + area(log(serBilir), IE_time = "IE_time"))
jm4a <- jm(cox, fm1, time_var = "year", functional_forms = ff4a)
jm4b <- jm(cox, fm1, time_var = "year", functional_forms = ff4b)
all(fixef(jm4a)[[1]] == fixef(jm4b)[[1]])


ff5a <- list("log(serBilir)" = ~ value(log(serBilir)) 
             + slope(log(serBilir), direction = "both")
             + area(log(serBilir)))
ff5b <- list("log(serBilir)" = ~ value(log(serBilir), IE_time = NULL) 
             + slope(log(serBilir), direction = "backward", IE_time = "IE_time")
             + area(log(serBilir), IE_time = "IE_time"))
jm5a <- jm(cox, fm1, time_var = "year", functional_forms = ff5a)
jm5b <- jm(cox, fm1, time_var = "year", functional_forms = ff5b)
all(fixef(jm5a)[[1]] == fixef(jm5b)[[1]])


ff6a <- list("log(serBilir)" = ~ value(log(serBilir)),
             "prothrombin" = ~ value(prothrombin))
ff6b <- list("log(serBilir)" = ~ value(log(serBilir), IE_time = "IE_time"),
             "prothrombin" = ~ value(prothrombin, IE_time = "IE_time2"))
jm6a <- jm(cox, list(fm1, fm2), time_var = "year", functional_forms = ff6a)
jm6b <- jm(cox, list(fm1, fm2), time_var = "year", functional_forms = ff6b)
all(fixef(jm6a)[[1]] == fixef(jm6b)[[1]])


ff7a <- list("log(serBilir)" = ~ Delta(log(serBilir)),
             "prothrombin" = ~ Delta(prothrombin))
ff7b <- list("log(serBilir)" = ~ Delta(log(serBilir), IE_time = "IE_time"),
             "prothrombin" = ~ Delta(prothrombin, IE_time = "IE_time2"))
jm7a <- jm(cox, list(fm1, fm2), time_var = "year", functional_forms = ff6a)
jm7b <- jm(cox, list(fm1, fm2), time_var = "year", functional_forms = ff6b)
all(fixef(jm7a)[[1]] == fixef(jm7b)[[1]])

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

fm1  <- lme(fixed = log(serBilir) ~ year + sex + age, 
            random =  ~ year | id, data = pbc2)
fm2 <- lme(prothrombin ~ year + sex, data = pbc2, random = ~ year | id)
pbc2.id$IE_time <- rep(Inf, 312)
cox <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)
ff7 <- list("log(serBilir)" = ~ value(log(serBilir), IE_time = "IE_time"),
            "prothrombin" = ~ value(prothrombin, IE_time = NULL))

jm7 <- jm(cox, list(fm1, fm2), time_var = "year", functional_forms = ff7)
all(fixef(jm1a)[[1]] == fixef(jm1b)[[1]])



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Surv_object <- cox 
# Mixed_objects <- list(fm1, fm2) 
# time_var <- "year"
# functional_forms <- ff7
# recurrent <- FALSE
# which_independent <- NULL
# base_hazard <- NULL
# data_Surv <- NULL
# id_var <- NULL
# priors <- NULL 
# control <- NULL
# source("R/help_functions.R")
# 
# 
# time <- st
# terms <- terms_FE_noResp
# data <- dataL # 'data' is the longitudinal dataset
# timeVar <- time_var 
# idVar <- idVar 
# idT <- idT # 'idT' is a character with the ID var name
# Fun_Forms <- collapsed_functional_forms
# Xbar <- Xbar
# eps <- eps 
# direction <- direction
# zero_ind <- zero_ind_X 
# time_window <- time_window
# IE_time <- IE_time

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

system("grep -R 'causal_effects' R/")

system("grep -R 'coefs' R/")
system("grep -R 'value' R/")
system("grep -R 'extract_attributes' R/")
system("grep -R 'design_matrices_functional_forms' R/")

system("grep -R 'desgn_matr' R/")
system("grep -R 'degn_matr_slp' R/")
system("grep -R 'degn_matr_acc' R/")
system("grep -R 'degn_matr_area' R/")

system("grep -R 'area' R/")
system("grep -R 'slope' R/")
system("grep -R 'velocity' R/")
system("grep -R 'acceleration' R/")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
rmarkdown::render("vignettes/JMbayes2.Rmd")
rmarkdown::render("vignettes/Baseline_Hazard.Rmd")
rmarkdown::render("vignettes/Causal_Effects.Rmd")
rmarkdown::render("vignettes/Competing_Risks.Rmd")
rmarkdown::render("vignettes/Dynamic_Predictions.Rmd")
rmarkdown::render("vignettes/Multi_State_Processes.Rmd")
rmarkdown::render("vignettes/Non_Gaussian_Mixed_Models.Rmd")
rmarkdown::render("vignettes/Recurring_Events.Rmd")
rmarkdown::render("vignettes/Time_Varying_Effects.Rmd")
rmarkdown::render("vignettes/Transformation_Functions.Rmd")
rmarkdown::render("vignettes/Super_Learning.Rmd")
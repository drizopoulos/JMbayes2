library('JMbayes2')
library("survival")
library("nlme")
library("GLMMadaptive")
library("coda")
library("lattice")
library("splines")
library("matrixStats")

fm1 <- lme(log(serBilir) ~ year * sex + I(year^2) + age + prothrombin,
           data = pbc2, random = ~ year | id)
fm2 <- lme(serChol ~ ns(year, 3) + sex + age, data = pbc2, random = ~ year | id,
           na.action = na.exclude)
fm3 <- mixed_model(hepatomegaly ~ sex + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())
fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())

#out1 <- log(pbc2$serBilir)
#length(out1)
#sum(is.na(out1))
#out2 <- pbc2$serChol
#length(out2)
#sum(is.na(out2))
#out3 <- pbc2$hepatomegaly
#length(out3)
#sum(is.na(out3))
#out4 <- pbc2$ascites
#length(out4)
#sum(is.na(out4))

#Mixed <- list(fm1, fm2, fm3, fm4)
Mixed <- list(fm1)
CoxFit <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

#save(CoxFit, Mixed, 
#     file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/reprod_fits_eta.RData')

obj <- jm(CoxFit, Mixed, time_var = "year")

JMbayes2::traceplot(obj, parm = 'betas')
JMbayes2:::ggtraceplot.jm

pbc2$prothrombin[pbc2$id == levels(pbc2$id)[1L]] <- NA
pbc2$prothrombin[pbc2$id == levels(pbc2$id)[2L]] <- NA
fm1 <- lme(log(serBilir) ~ year * (drug + sex) + I(year^2) + age + serChol,
           data = pbc2, random = ~ year + I(year^2)| id, na.action = na.exclude)
fm2 <- lme(prothrombin ~ ns(year, 2) + sex, data = pbc2,
           random = ~ year + I(year^2)| id,
           na.action = na.exclude, control = lmeControl(opt = "optim"))
Mixed <- list(fm1, fm2)
Cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)
obj <- jm(Cox, fm1, time_var = "year")

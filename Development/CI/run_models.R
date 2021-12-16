library("JMbayes2")
library("lattice")
source("./Development/CI/prepare_data.R")

#lmeFit <- lme(log(serBilir) ~ ns(year, k = 8, B = c(0, 14.5)) +
#                  IE + IE:ns(year - S, k = 4.5, B = c(0, 13.5)), data = pbc2,
#              random = ~ ns(year, k = 8, B = c(0, 14.5)) + IE +
#                  IE:ns(year - S, k = 4.5, B = c(0, 13.5)) | id,
#              control = lmeControl(opt = "optim"))

lmeFit <- lme(log(serBilir) ~ year + IE + IE:year, data = pbc2,
              random = ~ year + IE + IE:year | id,
              control = lmeControl(opt = "optim"))

CoxFit <- coxph(Surv(start, stop, event) ~ strata(CR) * (age + sex + IE),
                data = pbc2_CR)

fForms <- ~ value(log(serBilir)) + value(log(serBilir)):IE + value(log(serBilir)):CR

jointFit <- jm(CoxFit, lmeFit, time_var = "year", functional_forms = fForms,
               priors = list(Tau_alphas = list(diag(1), diag(1), 2000 * diag(1))),
               n_iter = 15000L, n_burnin = 5000L)

summary(jointFit)

lmeFit <- lme(log(serBilir) ~ year, data = pbc2, random = ~ year | id)

CoxFit <- coxph(Surv(years, status2) ~ age + sex, data = pbc2.id)

jointFit2 <- jm(CoxFit, lmeFit, time_var = "year")


save(list = c("jointFit", "jointFit2"),
     file = "./Development/CI/model.RData")




# the data frame that contains the combination of values to
# create the plot
for (tt in c(2, 3, 4, 5, 6, 7, 8, 9)) {
    t0 <- tt
    newDF <- with(pbc2, expand.grid(year = seq(0, 12, length.out = 55), S = t0))
    newDF$IE <- as.numeric(newDF$year >= newDF$S)

    # the effects plot
    print(xyplot(pred + low + upp ~ year,
                 data = effectPlotData(lmeFit, newDF, pbc2),
                 lty = c(1, 2, 2), col = c(2, 1, 1), lwd = 2, type = "l",
                 abline = list(v = tt, col = "blue", lty = 3, lwd = 3),
                 xlab = "Follow-up time (years)",
                 ylab = "log Serum Billirubin"))
}



Surv_object = CoxFit
Mixed_objects = lmeFit
time_var = 'year'
functional_forms = fForms
recurrent = FALSE
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL


lmeFit <- lme(log(serBilir) ~ ns(year, 2, B = c(0, 15)), data = pbc2,
              random = ~ ns(year, 2, B = c(0, 15)) | id)

CoxFit <- coxph(Surv(start, stop, event) ~ strata(CR) * (age + sex),
                data = pbc2_CR)

fForms <- ~ value(log(serBilir)) + value(log(serBilir)):CR

jointFit <- jm(CoxFit, lmeFit, time_var = "year", functional_forms = fForms,
               n_iter = 15000L, n_burnin = 5000L)

summary(jointFit)


newdataL <- pbc2[pbc2$id %in% c(81), ]
newdataL$status2 <- 0
newdataE <- pbc2_CR[pbc2_CR$id %in% c(81), ]
newdataE$event <- 0
newdataE <- newdataE[c(1, 3), ]
newdata <- list(newdataL = newdataL, newdataE = newdataE)

n <- nrow(newdataL)
for (i in seq_len(n)[-1]) {
   t0 <- newdataL$year[i]
   newdataL_i <- newdataL[newdataL$year <= t0, ]
   newdataL_i$years <- t0
   newdataE_i <- newdataE
   newdataE_i$stop <- t0 + 1e-06
   newdata_i <- list(newdataL = newdataL_i, newdataE = newdataE_i)
   ###
   predLong <- predict(jointFit, newdata = newdata_i,
                       times = seq(t0, 15, length.out = 51),
                       return_newdata = TRUE)
   predSurv <- predict(jointFit, newdata = newdata_i, process = "event",
                       return_newdata = TRUE)
   plot(predLong, predSurv, subject = 1)
}




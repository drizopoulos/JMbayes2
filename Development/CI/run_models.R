library("JMbayes2")
library("lattice")
source("./Development/CI/prepare_data.R")

lmeFit <- lme(log(serBilir) ~ ns(year, k = 8, B = c(0, 14.5)) +
                  IE + IE:ns(year - S, k = 4.5, B = c(0, 13.5)), data = pbc2,
              random = ~ ns(year, k = 8, B = c(0, 14.5)) + IE +
                  IE:ns(year - S, k = 4.5, B = c(0, 13.5)) | id,
              control = lmeControl(opt = "optim"))

lmeFit <- lme(log(serBilir) ~ year + IE + IE:year, data = pbc2,
              random = ~ year + IE + IE:year | id,
              control = lmeControl(opt = "optim"))

CoxFit <- coxph(Surv(start, stop, event) ~ strata(CR) * (age + sex), data = pbc2_CR)

fForms <- ~ value(log(serBilir)) + value(log(serBilir)):CR

jointFit <- jm(CoxFit, lmeFit, time_var = "year", functional_forms = fForms,
               priors = list(Tau_alphas = list(diag(1), 2000 * diag(1))),
               n_iter = 15000L, n_burnin = 5000L)

summary(jointFit)

save(list = "jointFit", file = "./Development/CI/model.RData")




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




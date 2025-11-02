library("JMbayes2")
CoxFit <- coxph(Surv(Time, death) ~ treat, data = prothros)

fm1 <- lme(pro ~ time, data = prothro, random = ~ time | id)
fm2 <- lme(pro ~ ns(time, 3), data = prothro,
           random = list(id = pdDiag(~ ns(time, 3))))
fm3 <- lme(pro ~ poly(time, 2), data = prothro,
           random = list(id = pdDiag(~ poly(time, 2))))

jointFit1 <- jm(CoxFit, fm1, time_var = "time")
jointFit2 <- jm(CoxFit, fm2, time_var = "time")
jointFit3 <- jm(CoxFit, fm3, time_var = "time")

Models <- list(jointFit1, jointFit2, jointFit3)
T0 <- 2.5
Data <- prothro[ave(prothro$time, prothro$id, FUN = max) > T0, ]
Data_after <- Data[Data$time > T0, ]
OptModel <- opt_model(Models, Data, T0, cores = 3L)
mises <- do.call('cbind', lapply(OptModel, function (x) {
    x[["std_MISEs_ave"]] + x[["std_MISEs_vario"]]
}))
best_model <- apply(mises, 1L, which.min)
Preds <- do.call('cbind', lapply(OptModel, function (x) x$Preds$newdata2$preds[[1L]]))
Obs <- Data_after$pro
id <- match(Data_after$id, unique(Data_after$id))

colMeans((Preds - Obs)^2)
mean((Preds[cbind(seq_along(id), best_model[id])] - Obs)^2)
weights <- t(apply(mises, 1L, function (x) exp(x) / sum(exp(x))))
mean((rowSums(weights[id, ] * Preds) - Obs)^2)

####

models = Models
t0 <- 3.5
newdata = Data
parallel = "snow"
cores = 1L

object = Models[[1]]

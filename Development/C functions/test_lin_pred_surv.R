library("JMbayes2")
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
fm2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)
fm3 <- mixed_model(ascites ~ year + sex, data = pbc2,
                   random = ~ year || id, family = binomial())

jointFit2 <- jm(CoxFit, list(fm1, fm2, fm3), time_var = "year",
                which_independent = cbind(1, 2))
summary(jointFit2)

X <- jointFit2$model_data$X_h
X[] <- lapply(X, function (x) unname(x[[1]]))
Z <- jointFit2$model_data$Z_h
Z[] <- lapply(Z, function (x) unname(x[[1]]))
betas <- list(jointFit2$statistics$Mean$betas1, jointFit2$statistics$Mean$betas2,
              jointFit2$statistics$Mean$betas3)
betas[] <- lapply(betas, unname)
b <- list(jointFit2$statistics$Mean$b[, 1:2], jointFit2$statistics$Mean$b[, 3:4],
          jointFit2$statistics$Mean$b[, 5:6])
id <- 0:311


tt1 <- linpred_surv(X, betas, Z, b, id)
tt2 <- linpred_surv2(X, betas, Z, b, id)

all.equal(tt1, tt2)

library("rbenchmark")

benchmark(
    Cpp1 = linpred_surv(X, betas, Z, b, id),
    Cpp2 = linpred_surv2(X, betas, Z, b, id),
    replications = 5000
)

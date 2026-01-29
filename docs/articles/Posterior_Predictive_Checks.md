# Posterior Predictive Checks

## Goodness-of-Fit Checks for Joint Models

To investigate the fit of a joint model to the training dataset \mathcal
D_n, we compare simulated longitudinal measurements \mathbf{y}\_i^{\sf
rep} and event times T_i^{\sf \*rep} from the fitted model with the
observed data \mathcal D_n. Ideally, we expect that realizations from
the fitted model are in close agreement with the observed data. This
framework supports multiple settings, including checks for existing
subjects, new subjects with only covariates, dynamic prediction at
intermediate follow-up times, and cross-validated assessment. For the
longitudinal component, goodness-of-fit is assessed through the mean,
variance, and correlation structure, while the survival component is
evaluated using empirical cumulative distributions and probability
integral transforms. The association between processes is examined using
time-dependent concordance statistics. A detailed presentation of this
framework is given in [Rizopoulos, Taylor and Kardys
(2026)](https://arxiv.org/abs/2601.18598).

## Posterior-Posterior Predictive Checks

``` r
CoxFit <- coxph(Surv(Time, death) ~ treat, data = prothros)

fm1 <- lme(pro ~ time * treat, data = prothro, random = ~ time | id)

jointFit1 <- jm(CoxFit, fm1, time_var = "time", save_random_effects = TRUE)
```

``` r
ppcheck(jointFit1, random_effects = "mcmc", type = "ecdf")
```

![](Posterior_Predictive_Checks_files/figure-html/Joint1-Long-eCDF-1.png)

``` r
ppcheck(jointFit1, random_effects = "mcmc", type = "average")
```

![](Posterior_Predictive_Checks_files/figure-html/Joint1-Long-mean-1.png)

``` r
ppcheck(jointFit1, random_effects = "mcmc", type = "variance")
```

![](Posterior_Predictive_Checks_files/figure-html/Joint1-Long-variance-1.png)

``` r
ppcheck(jointFit1, random_effects = "mcmc", type = "variogram")
```

![](Posterior_Predictive_Checks_files/figure-html/Joint1-Long-variogram-1.png)

``` r
FF1 <- function (t, betas, bi, data) {
    treat <- as.numeric(data$treat == "prednisone")
    X <- cbind(1, t, treat, t * treat)
    Z <- cbind(1, t)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}
```

``` r
ppcheck(jointFit1, random_effects = "mcmc", process = "event", Fforms_fun = FF1,
        type = "ecdf")
```

![](Posterior_Predictive_Checks_files/figure-html/Joint1-Surv-eCDF-1.png)

``` r
ppcheck(jointFit1, random_effects = "mcmc", process = "event", Fforms_fun = FF1,
        type = "surv-uniform")
```

![](Posterior_Predictive_Checks_files/figure-html/Joint1-Surv-ProbInt-1.png)

``` r
ppcheck(jointFit1, random_effects = "mcmc", process = "joint", Fforms_fun = FF1)
```

![](Posterior_Predictive_Checks_files/figure-html/Joint1-joint-1.png)

## Posterior-Prior Predictive Checks

``` r
fm2 <- lme(pro ~ ns(time, 3) * treat, data = prothro, 
           random = list(id = pdDiag(~ ns(time, 3))))

jointFit2 <- jm(CoxFit, fm2, time_var = "time")
```

``` r
FF2 <- function (t, betas, bi, data) {
    treat <- as.numeric(data$treat == "prednisone")
    NS <- ns(t, k = c(0.4928, 2.1547), B = c(0, 11.1078))
    X <- cbind(1, NS, treat, NS * treat)
    Z <- cbind(1, NS)
    eta <- c(X %*% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}
```

``` r
ppcheck(jointFit2, random_effects = "prior", type = "ecdf", Fforms_fun = FF2)
```

![](Posterior_Predictive_Checks_files/figure-html/Joint2-Long-1.png)

``` r
ppcheck(jointFit2, random_effects = "prior", type = "average", Fforms_fun = FF2)
```

![](Posterior_Predictive_Checks_files/figure-html/Joint2-Long-2.png)

``` r
ppcheck(jointFit2, random_effects = "prior", type = "variance", Fforms_fun = FF2)
```

![](Posterior_Predictive_Checks_files/figure-html/Joint2-Long-3.png)

``` r
ppcheck(jointFit2, random_effects = "prior", type = "variogram", Fforms_fun = FF2)
```

![](Posterior_Predictive_Checks_files/figure-html/Joint2-Long-4.png)

``` r
ppcheck(jointFit2, random_effects = "prior", process = "event", Fforms_fun = FF2,
        type = "ecdf")
```

![](Posterior_Predictive_Checks_files/figure-html/Joint2-Surv-1.png)

``` r
ppcheck(jointFit2, random_effects = "prior", process = "event", Fforms_fun = FF2,
        type = "surv-uniform")
```

![](Posterior_Predictive_Checks_files/figure-html/Joint2-Surv-2.png)

## Dynamic-Posterior-Posterior Predictive Checks

``` r
t0 <- 3
prothro_t0 <- prothro[prothro$Time > t0 & prothro$time <= t0, ]
prothro_t0$Time <- t0
prothro_t0$death <- 0
test_prothro <- prothro[prothro$Time > t0 & prothro$time > t0, ]
```

``` r
preds <- predict(jointFit2, newdata = prothro_t0, return_params_mcmc = TRUE)
```

``` r
ppcheck(jointFit2, random_effects = "mcmc", type = "average", 
        newdata = test_prothro, params_mcmc = preds$params_mcmc)
```

![](Posterior_Predictive_Checks_files/figure-html/Joint2-Long-dynamic-1.png)

## Cross-Validated Checks

``` r
CVdats <- create_folds(prothro, V = 5, id_var = "id")
```

``` r
fit_model <- function (data) {
    library("JMbayes2")
    data_id <- data[!duplicated(data$id), ]
    CoxFit <- coxph(Surv(Time, death) ~ treat, data = data_id)
    fm <- lme(pro ~ ns(time, 3) * treat, data = data, 
           random = list(id = pdDiag(~ ns(time, 3))))
    jm(CoxFit, fm, time_var = "time")
}

cl <- parallel::makeCluster(5L)
Models <- parallel::parLapply(cl, CVdats$training, fit_model)
parallel::stopCluster(cl)
```

``` r
ppcheck(Models, type = "variance", newdata = CVdats$testing,
        random_effects = "prior", Fforms_fun = FF2)
```

![](Posterior_Predictive_Checks_files/figure-html/CV_check-1.png)

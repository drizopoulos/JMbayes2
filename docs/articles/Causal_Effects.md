# Causal Effects

## Causal Effects from Joint Models

We will illustrate the calculation of causal effects from joint models
using the PBC dataset for the longitudinal outcome `serBilir` and the
composite event transplantation or death. We start by fitting a joint
model to the data. In the longitudinal submodel, we specify nonlinear
subject-specific trajectories using natural cubic splines. In the
fixed-effects part, we also include the treatment effect and its
interaction with time. In the survival submodel, we only include the
treatment effect.

``` r

pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")
lmeFit <- lme(log(serBilir) ~ ns(year, 3, B = c(0, 14.4)) * drug, 
                   data = pbc2, random = ~ ns(year, 3, B = c(0, 14.4)) | id,
                   control = lmeControl(opt = "optim"))
CoxFit <- coxph(Surv(years, status2) ~ drug, data = pbc2.id)
jmFit <- jm(CoxFit, lmeFit, time_var = "year")
summary(jmFit)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = CoxFit, Mixed_objects = lmeFit, time_var = "year")
#> 
#> Data Descriptives:
#> Number of groups: 312        Number of events: 169 (54.2%)
#> Number of observations:
#>   log(serBilir): 1945
#> 
#>                  DIC     WAIC      LPML
#> marginal    4227.440 4517.250 -2906.831
#> conditional 6550.666 6245.367 -3453.906
#> 
#> Random-effects covariance matrix:
#>                                                                       
#>                    StdDev   Corr                                      
#> (Intr)             0.9933 (Intr) n(,3,B=c(0,14.4))1 n(,3,B=c(0,14.4))2
#> n(,3,B=c(0,14.4))1 1.5480 0.2289                                      
#> n(,3,B=c(0,14.4))2 1.6884 0.4480 0.7257                               
#> n(,3,B=c(0,14.4))3 1.8900 0.4437 0.1676             0.6655            
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%      P   Rhat
#> drugD-penicil        -0.0175 0.2040 -0.4034 0.3710 0.9456 1.0133
#> value(log(serBilir))  1.2912 0.0903  1.1192 1.4693 0.0000 1.0598
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                        Mean  StDev    2.5%  97.5%      P   Rhat
#> (Intercept)          0.5877 0.0821  0.4261 0.7491 0.0000 1.0019
#> ns(,3,B=c(0,14.4))1  1.1569 0.1751  0.8200 1.5037 0.0000 1.0531
#> ns(,3,B=c(0,14.4))2  2.1556 0.2216  1.7403 2.6023 0.0000 1.2445
#> ns(,3,B=c(0,14.4))3  2.3385 0.3193  1.7538 2.9917 0.0000 1.3621
#> drugD-penicil       -0.1078 0.1162 -0.3365 0.1188 0.3502 1.0009
#> n(,3,B=c(0,14.4))1:  0.1781 0.2388 -0.2942 0.6511 0.4598 1.0039
#> n(,3,B=c(0,14.4))2: -0.4587 0.2848 -1.0075 0.1117 0.1078 1.0621
#> n(,3,B=c(0,14.4))3: -0.6746 0.4060 -1.4589 0.1376 0.0978 1.1464
#> sigma                0.2886 0.0063  0.2770 0.3012 0.0000 1.0063
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 16 sec
```

The coefficient for `drugD-penicil` for the survival outcome in the
output produced by the
[`summary()`](https://rdrr.io/r/base/summary.html) method denotes the
residual/direct effect of treatment on the risk of the composite event.
It does not include the effect of treatment that follows via the serum
bilirubin pathway.

We will illustrate the calculation of causal risk differences for the
group of patients that have the same distribution of serum bilirubin
values as Patient 2:

``` r

xyplot(log(serBilir) ~ year, data = pbc2, subset = id == 2, type = "b",
       xlab = "Follow-up time (years)", ylab = "log{serum bilirubin (mg/dL)}",
       main = "Patient 2")
```

![](Causal_Effects_files/figure-html/unnamed-chunk-1-1.png)

We calculate the risk difference for the composite event between the
active treatment D-penicillamine and placebo at the horizon time
`t_horiz = 6` using the longitudinal data up to year `t0 = 4`. To
achieve this, we create a dataset with this patient’s data. This patient
received the active treatment D-penicillamine; hence, we also create a
version of her data with the `drug` variable set to `placebo`:

``` r

t0 <- 4
t_horiz <- 6
dataP2_Dpenici <- pbc2[pbc2$id == 2 & pbc2$year <= t0, ]
dataP2_Dpenici$years <- t0
dataP2_Dpenici$status2 <- 0

dataP2_placebo <- dataP2_Dpenici
dataP2_placebo$drug <- factor("placebo", levels = levels(pbc2$drug))
```

Note that in the `dataP2_placebo` dataset, we need to specify that
`drug` is a factor with two levels. We also specify that the last time
point we know the patient was still event-free was `t0`.

We estimate the cumulative risk for the composite event at `t_horiz`
under the active treatment arm using the
[`predict()`](https://rdrr.io/r/stats/predict.html) method:

``` r

Pr1 <- predict(jmFit, newdata = dataP2_Dpenici, process = "event", 
               times = t_horiz, return_mcmc = TRUE)
```

We have set the argument `return_mcmc` to `TRUE` to enable the
calculation of a credible interval that accounts for the MCMC
uncertainty. We produce the same estimate under the placebo arm:

``` r

Pr0 <- predict(jmFit, newdata = dataP2_placebo, process = "event", 
               times = t_horiz, return_mcmc = TRUE)
```

The estimated risk difference and its 95% credible interval are
calculated by the corresponding elements of the `Pr1` and `Pr0` objects,
i.e.,

``` r

# estimate 
Pr1$pred[2L] - Pr0$pred[2L]
#> [1] 0.0001722797

# MCMC variability
quantile(Pr1$mcmc[2L, ] - Pr0$mcmc[2L, ], probs = c(0.025, 0.975))
#>       2.5%      97.5% 
#> -0.2249485  0.2152505
```

### Time-varying treatments

An extended example with time-varying treatments / intermediate events
that showcases a calculation of the variance of the causal effects that
includes the sampling variability is available
[here](https://github.com/drizopoulos/JMbayes2/blob/master/Development/CI/causal_effects.R).

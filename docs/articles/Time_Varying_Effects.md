# Time Varying Effects

## Non Proportional Hazards

The basic definition of the joint model assumes the coefficients that
quantify the association between the versions of the longitudinal
outcomes and the hazard of the event are time-constant (i.e., the
proportional hazards assumption). We can relax this assumption by
specifying time-varying coefficients via the `functional_forms` argument
of function
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md).

We will illustrate this capability using the PBC dataset. We start by
fitting a Cox model for the composite event transplantation or death,
including sex as a baseline covariate:

``` r

pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
```

We aim to assess the strength of the association between the risk of the
composite event and the serum bilirubin level. We will describe the
patient-specific profiles over time for this biomarker using a linear
mixed-effects model, where we include an intercept in both the fixed and
random effects, as well as the linear and quadratic time effects. In the
fixed effects, we also include the interaction of the time effect and
sex. The syntax to fit this model with
[`lme()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
is:

``` r

fm <- lme(log(serBilir) ~ poly(year, 2) * sex, data = pbc2, 
          random = ~ poly(year, 2) | id, control = lmeControl(opt = 'optim'))
```

The default call to
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md) adds
the subject-specific linear predictor of the mixed model as a
time-varying covariate in the survival relative risk model:

``` r

jointFit1 <- jm(CoxFit, fm, time_var = "year")
summary(jointFit1)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = CoxFit, Mixed_objects = fm, time_var = "year")
#> 
#> Data Descriptives:
#> Number of groups: 312        Number of events: 169 (54.2%)
#> Number of observations:
#>   log(serBilir): 1945
#> 
#>                  DIC     WAIC      LPML
#> marginal    4344.764 6074.567 -3295.282
#> conditional 8720.740 8458.352 -4517.665
#> 
#> Random-effects covariance matrix:
#>                               
#>        StdDev    Corr         
#> (Intr) 1.3098  (Intr)  p(,2)1 
#> p(,2)1 21.6589 0.6792         
#> p(,2)2 12.2916 -0.2326 -0.1198
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale            -0.1625 0.2696 -0.6530 0.4066 0.5307 1.0052
#> value(log(serBilir))  1.2986 0.0963  1.1241 1.4966 0.0000 1.0733
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev     2.5%   97.5%      P   Rhat
#> (Intercept)     1.4871 0.2261   1.0356  1.9304 0.0000 1.0019
#> poly(year, 2)1 29.3995 5.1220  19.4766 39.6132 0.0000 1.0290
#> poly(year, 2)2 -4.8924 3.1649 -11.4537  1.0405 0.1093 1.0134
#> sexfemale      -0.4702 0.2398  -0.9366  0.0037 0.0513 1.0002
#> p(,2)1         -4.9553 5.3984 -15.6122  5.7616 0.3524 1.0317
#> p(,2)2          6.3926 3.3626   0.0091 13.2588 0.0489 1.0196
#> sigma           0.3028 0.0062   0.2910  0.3148 0.0000 1.0174
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 15 sec
```

To specify that the association of serum bilirubin may change over time,
we include an interaction of this time-varying covariate with a natural
cubic spline of time using function
[`ns()`](https://rdrr.io/r/splines/ns.html) from the **splines**
package. **Important Note:** For this to work correctly, we need to
explicitly specify the internal and boundary knots for the B-splines
basis, i.e., in the following example, we set the internal knots at 3,
6, and 9 years, and the boundary knots at 0 and 14.5 years:

``` r

form_splines <- ~ value(log(serBilir)) * ns(year, k = c(3, 6, 9), B = c(0, 14.5))
jointFit2 <- update(jointFit1, functional_forms = form_splines, 
                    n_iter = 7000L, n_burnin = 1500L)
summary(jointFit2)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = CoxFit, Mixed_objects = fm, time_var = "year", 
#>     functional_forms = form_splines, n_iter = 7000L, n_burnin = 1500L)
#> 
#> Data Descriptives:
#> Number of groups: 312        Number of events: 169 (54.2%)
#> Number of observations:
#>   log(serBilir): 1945
#> 
#>                  DIC     WAIC      LPML
#> marginal    4365.120 6286.215 -3598.082
#> conditional 8718.259 8454.858 -4546.904
#> 
#> Random-effects covariance matrix:
#>                               
#>        StdDev    Corr         
#> (Intr) 1.3109  (Intr)  p(,2)1 
#> p(,2)1 21.8692 0.6847         
#> p(,2)2 12.0973 -0.2244 -0.1201
#> 
#> Survival outcome:
#>                                                                   Mean  StDev
#> sexfemale                                                      -0.1753 0.2758
#> value(log(serBilir))                                            1.3700 0.2432
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))1 -0.2646 0.3003
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))2  0.0100 0.3451
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))3 -0.5594 0.7105
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))4 -1.1220 0.8765
#>                                                                   2.5%  97.5%
#> sexfemale                                                      -0.6922 0.3753
#> value(log(serBilir))                                            0.9336 1.8889
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))1 -0.8700 0.2980
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))2 -0.6837 0.6898
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))3 -1.9819 0.7919
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))4 -2.8159 0.5355
#>                                                                     P   Rhat
#> sexfemale                                                      0.5179 1.0004
#> value(log(serBilir))                                           0.0000 1.0291
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))1 0.3782 1.0330
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))2 0.9822 1.0257
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))3 0.4570 1.0635
#> value(log(serBilir)):ns(year, k = c(3, 6, 9), B = c(0, 14.5))4 0.1967 1.0727
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev     2.5%   97.5%      P   Rhat
#> (Intercept)     1.4937 0.2250   1.0545  1.9375 0.0000 1.0000
#> poly(year, 2)1 29.7802 5.2619  19.7385 40.1633 0.0000 1.0033
#> poly(year, 2)2 -4.5358 3.1729 -10.8049  1.6951 0.1487 1.0059
#> sexfemale      -0.4734 0.2397  -0.9417 -0.0001 0.0499 1.0002
#> p(,2)1         -5.1671 5.4772 -15.7512  5.5301 0.3468 1.0171
#> p(,2)2          6.0727 3.2851  -0.3268 12.5543 0.0606 1.0356
#> sigma           0.3031 0.0061   0.2916  0.3155 0.0000 1.0038
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 7000 
#> burn-in per chain: 1500 
#> thinning: 1 
#> time: 28 sec
```

The spline coefficients do not have a straightforward interpretation.
We, therefore, visualize the time-varying association of log serum
bilirubin with the hazard of the composite event using the following
piece of code:

``` r

x_times <- seq(0.001, 12, length = 501)
X <- cbind(1, ns(x_times, knots = c(3, 6, 9), B = c(0, 14.5)))
mcmc_alphas <- do.call('rbind', jointFit2$mcmc$alphas)
log_hr <- X %*% t(mcmc_alphas)
log_hr_mean <- rowMeans(log_hr)
log_hr_low <- apply(log_hr, 1, quantile, probs = 0.025)
log_hr_upp <- apply(log_hr, 1, quantile, probs = 0.975)

matplot(x_times, cbind(exp(log_hr_mean), exp(log_hr_low), exp(log_hr_upp)), 
        type = "l", col = c("red", "black", "black"), lty = c(1, 2, 2), lwd = 2,
        xlab = "Follow-up Time (years)", ylab = "Hazard Ratio log serum Bilirubin",
        ylim = c(0.5, 6.4))
abline(h = exp(coef(jointFit1)$association), lty = 2, col = "red")
abline(h = 1, lty = 2)
legend("topright", c("time-varying coefficient", "proportional hazards"),
       lty = c(1, 2), lwd = c(2, 1), col = "red", bty = "n")
```

![](Time_Varying_Effects_files/figure-html/unnamed-chunk-5-1.png)

We observe that the 95% credible interval for the time-varying
coefficient includes the horizontal line corresponding to proportional
hazards. This is also confirmed by comparing the two models:

``` r

compare_jm(jointFit1, jointFit2)
#> 
#>                 DIC     WAIC      LPML
#>  jointFit1 4344.764 6074.567 -3295.282
#>  jointFit2 4365.120 6286.215 -3598.082
#> 
#> The criteria are calculated on the basis of the marginal log-likelihood.
```

The WAIC and LPML indicate that `jointFit1` is a better model than
`jointFit2`. The DIC has the same magnitude for both models.

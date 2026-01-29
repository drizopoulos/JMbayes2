# Univariate and Multivariate Joint Models

## Fitting Joint Models with JMbayes2

### Univariate

The function that fits joint models in **JMbayes2** is called
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md). It has
three required arguments, `Surv_object` a Cox model fitted by
[`coxph()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
or an Accelerated Failure time model fitted by
[`survreg()`](https://rdrr.io/pkg/survival/man/survreg.html),
`Mixed_objects` a single or a list of mixed models fitted either by the
[`lme()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
or
[`mixed_model()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
functions, and `time_var` a character string indicating the name of the
time variable in the specification of the mixed-effects models. We will
illustrate the basic use of the package in the PBC dataset. We start by
fitting a Cox model for the composite event transplantation or death,
including sex as a baseline covariate:

``` r
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
```

We aim to assess the strength of the association between the risk of the
composite event and the serum bilirubin levels collected during
follow-up. We will describe the patient-specific profiles over time for
this biomarker using a linear mixed model, with fixed effects, time,
sex, and their interaction, and as random effects, random intercepts,
and random slopes. The syntax to fit this model with
[`lme()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
is:

``` r
fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
```

The joint model that links the survival and longitudinal submodels is
fitted with the following call to the
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md)
function:

``` r
jointFit1 <- jm(CoxFit, fm1, time_var = "year")
summary(jointFit1)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = CoxFit, Mixed_objects = fm1, time_var = "year")
#> 
#> Data Descriptives:
#> Number of Groups: 312        Number of events: 169 (54.2%)
#> Number of Observations:
#>   log(serBilir): 1945
#> 
#>                  DIC     WAIC      LPML
#> marginal    4361.435 5361.220 -3356.241
#> conditional 3536.629 3355.317 -1907.678
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 1.0028 (Intr)
#> year   0.1829 0.3994
#> 
#> Survival Outcome:
#>                         Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale            -0.1581 0.2717 -0.6499 0.3848 0.5544 1.0015
#> value(log(serBilir))  1.2433 0.0847  1.0776 1.4140 0.0000 1.0183
#> 
#> Longitudinal Outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)     0.7239 0.1720  0.3821  1.0600 0.0000 0.9997
#> year            0.2668 0.0381  0.1929  0.3444 0.0000 1.0024
#> sexfemale      -0.2639 0.1823 -0.6192  0.0882 0.1511 0.9999
#> year:sexfemale -0.0886 0.0404 -0.1681 -0.0093 0.0247 1.0028
#> sigma           0.3465 0.0065  0.3342  0.3596 0.0000 1.0101
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 17 sec
```

The output of the [`summary()`](https://rdrr.io/r/base/summary.html)
method provides some descriptive statistics of the sample at hand,
followed by some fit statistics based on the marginal (random effects
are integrated out using the Laplace approximation) and conditional on
the random effects log-likelihood functions, followed by the estimated
variance-covariance matrix for the random effects, followed by the
estimates for the survival submodel, followed by the estimates for the
longitudinal submodel(s), and finally some information for the MCMC
fitting algorithm.

By default,
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md) adds
the subject-specific linear predictor of the mixed model as a
time-varying covariate in the survival relative risk model. In the
output, this is named as `value(log(serBilir))` to denote that, by
default, the current value functional form is used. That is, we assume
that the instantaneous risk of an event at a specific time t is
associated with the value of the linear predictor of the longitudinal
outcome at the same time point t.

Standard MCMC diagnostics are available to evaluate convergence. For
example, the traceplot for the association coefficient
`value(log(serBilir))` is produced with the following syntax:

``` r
ggtraceplot(jointFit1, "alphas")
```

![](JMbayes2_files/figure-html/unnamed-chunk-4-1.png)

and the density plot with the call:

``` r
ggdensityplot(jointFit1, "alphas")
```

![](JMbayes2_files/figure-html/unnamed-chunk-5-1.png)

#### Notes

- The ordering of the subjects in the datasets used to fit the mixed and
  Cox regression models needs to be the same.

- The units of the time variables in the mixed and Cox models need to be
  the same.

### Multivariate

To fit a joint model with multiple longitudinal outcomes, we provide a
list of mixed models as the second argument of
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md). In the
following example, we extend the joint model we fitted above by
including the prothrombin time and the log odds of the presence or
absence of ascites as time-varying covariates in the relative risk model
for the composite event. Ascites is a dichotomous outcome, and
therefore, we fit a mixed-effects logistic regression model for it using
the
[`mixed_model()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
function from the **GLMMadaptive** package. The use of `||` in the
`random` argument of
[`mixed_model()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
specifies that the random intercepts and random slopes are assumed
uncorrelated. In addition, the argument `which_independent` can be used
to determine which longitudinal outcomes are to be assumed independent;
here, as an illustration, we specify that the first (i.e., serum
bilirubin) and second (i.e., prothrombin time) longitudinal outcomes are
independent. To assume that all longitudinal outcomes are independent,
we can use `jm(..., which_independent = "all")`. Because this joint
model is more complex, we increase the number of MCMC iterations, the
number of burn-in iterations, and the thinning per chain using the
corresponding control arguments:

``` r
fm2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)
fm3 <- mixed_model(ascites ~ year + sex, data = pbc2,
                   random = ~ year || id, family = binomial())

jointFit2 <- jm(CoxFit, list(fm1, fm2, fm3), time_var = "year",
                which_independent = cbind(1, 2),
                n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)
summary(jointFit2)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = CoxFit, Mixed_objects = list(fm1, 
#>     fm2, fm3), time_var = "year", which_independent = cbind(1, 
#>     2), n_iter = 12000L, n_burnin = 2000L, n_thin = 5L)
#> 
#> Data Descriptives:
#> Number of Groups: 312        Number of events: 169 (54.2%)
#> Number of Observations:
#>   log(serBilir): 1945
#>   prothrombin: 1945
#>   ascites: 1885
#> 
#>                  DIC     WAIC      LPML
#> marginal    11655.22 16089.26 -8730.453
#> conditional 12879.91 12590.61 -6812.813
#> 
#> Random-effects covariance matrix:
#>                                                   
#>        StdDev   Corr                              
#> (Intr) 1.0022 (Intr)   year (Intr)    year  (Intr)
#> year   0.1866 0.4490                              
#> (Intr) 0.7625                                     
#> year   0.3241               -0.0122               
#> (Intr) 2.7049 0.5177 0.4745 0.3283  -0.0298       
#> year   0.4613 0.4057 0.6660 -0.0592 0.3448        
#> 
#> Survival Outcome:
#>                         Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale            -0.6621 0.3613 -1.3655 0.0338 0.0607 1.0140
#> value(log(serBilir))  0.4863 0.1786  0.1096 0.8212 0.0147 1.0545
#> value(prothrombin)   -0.0583 0.1244 -0.3296 0.1735 0.6293 1.0612
#> value(ascites)        0.6227 0.1460  0.3708 0.9518 0.0000 1.0703
#> 
#> Longitudinal Outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%     P   Rhat
#> (Intercept)     0.6926 0.1691  0.3584  1.0311 0.000 1.0003
#> year            0.2694 0.0349  0.2005  0.3383 0.000 1.0004
#> sexfemale      -0.2357 0.1795 -0.5953  0.1183 0.190 1.0002
#> year:sexfemale -0.0800 0.0362 -0.1508 -0.0097 0.024 1.0022
#> sigma           0.3480 0.0068  0.3347  0.3617 0.000 1.0047
#> 
#> Longitudinal Outcome: prothrombin (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    10.9863 0.1728 10.6532 11.3254 0.0000 1.0033
#> year            0.2081 0.0774  0.0592  0.3599 0.0070 1.0065
#> sexfemale      -0.4422 0.1831 -0.8040 -0.0912 0.0100 1.0051
#> year:sexfemale  0.0470 0.0809 -0.1130  0.2029 0.5577 1.0088
#> sigma           1.0569 0.0202  1.0185  1.0975 0.0000 1.0004
#> 
#> Longitudinal Outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept) -4.4912 0.6735 -5.9197 -3.2356 0.0000 1.0121
#> year         0.6393 0.0684  0.5128  0.7854 0.0000 1.0651
#> sexfemale   -0.5556 0.6565 -1.8222  0.7787 0.3913 1.0021
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 12000 
#> burn-in per chain: 2000 
#> thinning: 5 
#> time: 1.9 min
```

The survival submodel output now contains the estimated coefficients for
`value(prothrombin)` and `value(ascites)`, as well as parameter
estimates for all three longitudinal submodels.

### Functional forms

As mentioned above, the default call to
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md)
includes the subject-specific linear predictors of the mixed-effects
models as time-varying covariates in the relative risk model. However,
this is just one of the many possibilities for linking longitudinal and
survival outcomes. The argument `functional_forms` of
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md)
provides additional options. Based on previous experience, two extra
functional forms are provided: the time-varying slope and the
time-varying *normalized* area/cumulative effect. The time-varying slope
is the first-order derivative of the subject-specific linear predictor
of the mixed-effect model with respect to the (follow-up) time variable.
The time-varying *normalized* area/cumulative effect is the integral of
the subject-specific linear predictor of the mixed-effect model from
zero to the current (follow-up) time t divided by t. The integral is the
area under the subject-specific longitudinal profile; by dividing the
integral by t, we obtain the average of the subject-specific
longitudinal profile over the corresponding period (0, t).

To illustrate how the `functional_forms` argument can be used to specify
these functional forms, we update the joint model `jointFit2` by
including the time-varying slope of log serum bilirubin instead of the
value and also the interaction of this slope with sex and for
prothrombin we include the normalized cumulative effect. For ascites, we
keep the current value functional form. The corresponding syntax to fit
this model is:

``` r
fForms <- list(
  "log(serBilir)" = ~ slope(log(serBilir)) + slope(log(serBilir)):sex,
  "prothrombin"   = ~ JMbayes2::area(prothrombin)
)

jointFit3 <- update(jointFit2, functional_forms = fForms)
summary(jointFit3)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = CoxFit, Mixed_objects = list(fm1, 
#>     fm2, fm3), time_var = "year", functional_forms = fForms, 
#>     which_independent = cbind(1, 2), n_iter = 12000L, n_burnin = 2000L, ...
#> 
#> Data Descriptives:
#> Number of Groups: 312        Number of events: 169 (54.2%)
#> Number of Observations:
#>   log(serBilir): 1945
#>   prothrombin: 1945
#>   ascites: 1885
#> 
#>                  DIC     WAIC      LPML
#> marginal    11692.23 12995.55 -7306.362
#> conditional 12656.91 12381.64 -6669.114
#> 
#> Random-effects covariance matrix:
#>                                                   
#>        StdDev   Corr                              
#> (Intr) 0.9989 (Intr)   year (Intr)    year  (Intr)
#> year   0.1853 0.4548                              
#> (Intr) 0.7500                                     
#> year   0.3232               -0.0047               
#> (Intr) 2.5559 0.5529 0.4692 0.3487  -0.0775       
#> year   0.4361 0.4303 0.6690 -0.0720 0.3773        
#> 
#> Survival Outcome:
#>                                   Mean  StDev     2.5%  97.5%      P   Rhat
#> sexfemale                       0.3595 0.9644  -1.3670 2.3829 0.7557 1.0662
#> slope(log(serBilir))            4.3938 2.5197  -0.0624 9.5851 0.0577 1.1294
#> slope(log(serBilir)):sexfemale -4.3880 3.0045 -11.1525 0.5844 0.1020 1.1396
#> JMbayes2::area(prothrombin)    -0.4097 0.3097  -0.9998 0.1627 0.1957 1.3476
#> value(ascites)                  1.0639 0.2373   0.6257 1.5698 0.0000 1.3015
#> 
#> Longitudinal Outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)     0.6684 0.1660  0.3397  0.9945 0.0003 1.0079
#> year            0.2658 0.0336  0.2025  0.3346 0.0000 1.0005
#> sexfemale      -0.2049 0.1768 -0.5497  0.1392 0.2433 1.0087
#> year:sexfemale -0.0745 0.0348 -0.1450 -0.0079 0.0270 1.0008
#> sigma           0.3483 0.0066  0.3352  0.3615 0.0000 1.0046
#> 
#> Longitudinal Outcome: prothrombin (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    10.9993 0.1684 10.6565 11.3294 0.0000 1.0003
#> year            0.1839 0.0764  0.0325  0.3362 0.0140 1.0056
#> sexfemale      -0.4574 0.1786 -0.7983 -0.1039 0.0140 1.0006
#> year:sexfemale  0.0702 0.0806 -0.0884  0.2279 0.3857 1.0058
#> sigma           1.0591 0.0203  1.0197  1.0994 0.0000 1.0152
#> 
#> Longitudinal Outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept) -4.4105 0.6235 -5.7030 -3.2274 0.0000 1.0378
#> year         0.6304 0.0777  0.4830  0.7938 0.0000 1.2194
#> sexfemale   -0.4225 0.6122 -1.6453  0.7616 0.4833 1.0057
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 12000 
#> burn-in per chain: 2000 
#> thinning: 5 
#> time: 2 min
```

As seen above, the `functional_forms` argument is a named list with
elements corresponding to the longitudinal outcomes. If a longitudinal
outcome is not specified in this list, then the default value functional
form is used for that outcome. Each element of the list should be a
one-sided R formula in which the functions
[`value()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md),
[`slope()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md), and
[`area()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md) can
be used. Interaction terms between the functional forms and other
(baseline) covariates are also allowed.

### Penalized Coefficients using Shrinkage Priors

When multiple longitudinal outcomes are considered with possibly
different functional forms per outcome, we require to fit a relative
risk model containing several terms. Moreover, it is often of scientific
interest to select which terms/functional forms per longitudinal outcome
are more strongly associated with the risk of the event of interest. To
facilitate this selection,
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md) allows
penalizing the regression coefficients using shrinkage priors. As an
example, we refit `jointFit3` by assuming a Horseshoe prior for the
`alphas` coefficients (i.e., the coefficients of the longitudinal
outcomes in the relative risk model):

``` r
jointFit4 <- update(jointFit3, priors = list("penalty_alphas" = "horseshoe"))
cbind("un-penalized" = unlist(coef(jointFit3)), 
      "penalized" = unlist(coef(jointFit4)))
#>                                            un-penalized  penalized
#> gammas.Mean                                   0.3594950 -0.5053443
#> association.slope(log(serBilir))              4.3938467  2.0175752
#> association.slope(log(serBilir)):sexfemale   -4.3879624 -0.8996426
#> association.JMbayes2::area(prothrombin)      -0.4096893 -0.1321933
#> association.value(ascites)                    1.0639313  0.8572684
```

Apart from the Horseshoe prior, the ridge prior is also provided.

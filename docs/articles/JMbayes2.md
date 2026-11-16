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
#> Number of groups: 312        Number of events: 169 (54.2%)
#> Number of observations:
#>   log(serBilir): 1945
#> 
#>                  DIC     WAIC      LPML
#> marginal    4366.988 5367.024 -2995.202
#> conditional 3542.137 3360.334 -1928.831
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 1.0028 (Intr)
#> year   0.1844 0.4087
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%    P   Rhat
#> sexfemale            -0.1628 0.2721 -0.6623 0.4168 0.53 1.0014
#> value(log(serBilir))  1.2536 0.0877  1.0859 1.4311 0.00 1.0250
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)     0.7216 0.1735  0.3749  1.0577 0.0002 1.0002
#> year            0.2668 0.0384  0.1917  0.3437 0.0000 1.0010
#> sexfemale      -0.2628 0.1849 -0.6235  0.1074 0.1553 1.0001
#> year:sexfemale -0.0877 0.0407 -0.1683 -0.0083 0.0320 1.0010
#> sigma           0.3470 0.0067  0.3341  0.3602 0.0000 1.0164
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 14 sec
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
#> Number of groups: 312        Number of events: 169 (54.2%)
#> Number of observations:
#>   log(serBilir): 1945
#>   prothrombin: 1945
#>   ascites: 1885
#> 
#>                  DIC     WAIC      LPML
#> marginal    11643.66 13683.69 -7469.797
#> conditional 12856.57 12597.17 -6801.781
#> 
#> Random-effects covariance matrix:
#>                                                   
#>        StdDev   Corr                              
#> (Intr) 1.0004 (Intr)   year (Intr)    year  (Intr)
#> year   0.1860 0.4422                              
#> (Intr) 0.7586                                     
#> year   0.3244               -0.0073               
#> (Intr) 2.7106 0.5163 0.4763 0.3309  -0.0260       
#> year   0.4576 0.3795 0.6547 -0.0553 0.3626        
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale            -0.6616 0.3728 -1.3902 0.0620 0.0730 1.0243
#> value(log(serBilir))  0.5221 0.1722  0.1772 0.8526 0.0027 1.1023
#> value(prothrombin)   -0.0600 0.1265 -0.3298 0.1626 0.6800 1.0403
#> value(ascites)        0.5971 0.1515  0.3530 0.9412 0.0000 1.2396
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%     P   Rhat
#> (Intercept)     0.7014 0.1673  0.3734  1.0318 0.000 1.0050
#> year            0.2702 0.0358  0.2005  0.3416 0.000 1.0007
#> sexfemale      -0.2434 0.1777 -0.5995  0.1023 0.175 1.0061
#> year:sexfemale -0.0811 0.0372 -0.1531 -0.0063 0.032 1.0031
#> sigma           0.3479 0.0068  0.3348  0.3612 0.000 1.0035
#> 
#> Longitudinal outcome: prothrombin (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    10.9817 0.1667 10.6529 11.3003 0.0000 1.0027
#> year            0.2113 0.0754  0.0663  0.3636 0.0077 1.0057
#> sexfemale      -0.4395 0.1776 -0.7882 -0.0892 0.0127 1.0036
#> year:sexfemale  0.0445 0.0786 -0.1110  0.1988 0.5670 1.0066
#> sigma           1.0581 0.0205  1.0183  1.0987 0.0000 1.0000
#> 
#> Longitudinal outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept) -4.4747 0.6986 -5.9347 -3.1790 0.0000 1.0810
#> year         0.6297 0.0765  0.4906  0.7902 0.0000 1.2771
#> sexfemale   -0.5556 0.6575 -1.8060  0.7473 0.3963 1.0086
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 12000 
#> burn-in per chain: 2000 
#> thinning: 5 
#> time: 1.2 min
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
#> Number of groups: 312        Number of events: 169 (54.2%)
#> Number of observations:
#>   log(serBilir): 1945
#>   prothrombin: 1945
#>   ascites: 1885
#> 
#>                  DIC     WAIC      LPML
#> marginal    11678.65 12627.22 -6967.448
#> conditional 12746.22 12468.45 -6732.847
#> 
#> Random-effects covariance matrix:
#>                                                    
#>        StdDev   Corr                               
#> (Intr) 0.9957 (Intr)   year (Intr)    year  (Intr) 
#> year   0.1858 0.4586                               
#> (Intr) 0.7557                                      
#> year   0.3239               -0.0096                
#> (Intr) 2.6274 0.5337 0.4527 0.3538  -0.0693        
#> year   0.4626 0.4526 0.6744 -0.0688 0.3579  -0.0000
#> 
#> Survival outcome:
#>                                   Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale                       0.2129 0.8797 -1.3963 2.0428 0.8333 1.0258
#> slope(log(serBilir))            4.3599 2.4154 -0.4874 8.8593 0.0780 1.0467
#> slope(log(serBilir)):sexfemale -3.8157 2.7106 -9.5174 1.1206 0.1543 1.0729
#> JMbayes2::area(prothrombin)    -0.3629 0.3291 -1.0444 0.1903 0.2653 1.3639
#> value(ascites)                  0.9871 0.2506  0.5929 1.5823 0.0000 1.5166
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)     0.6594 0.1661  0.3336  0.9778 0.0000 1.0007
#> year            0.2637 0.0340  0.1988  0.3310 0.0000 1.0068
#> sexfemale      -0.1971 0.1756 -0.5362  0.1488 0.2623 1.0003
#> year:sexfemale -0.0718 0.0349 -0.1405 -0.0037 0.0407 1.0038
#> sigma           0.3481 0.0066  0.3352  0.3613 0.0000 1.0017
#> 
#> Longitudinal outcome: prothrombin (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    10.9939 0.1709 10.6604 11.3302 0.0000 1.0020
#> year            0.1896 0.0771  0.0383  0.3416 0.0130 1.0108
#> sexfemale      -0.4512 0.1792 -0.8083 -0.1001 0.0123 1.0012
#> year:sexfemale  0.0641 0.0802 -0.0965  0.2224 0.4200 1.0127
#> sigma           1.0586 0.0201  1.0199  1.0991 0.0000 1.0065
#> 
#> Longitudinal outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept) -4.5441 0.6871 -5.9249 -3.2525 0.0000 1.0596
#> year         0.6486 0.0788  0.4999  0.7970 0.0000 1.2031
#> sexfemale   -0.3832 0.6311 -1.6099  0.8633 0.5307 1.0011
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 12000 
#> burn-in per chain: 2000 
#> thinning: 5 
#> time: 1.2 min
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
#> gammas.Mean                                   0.2129284 -0.5296200
#> association.slope(log(serBilir))              4.3599177  1.9665255
#> association.slope(log(serBilir)):sexfemale   -3.8157407 -0.7993484
#> association.JMbayes2::area(prothrombin)      -0.3628543 -0.1529725
#> association.value(ascites)                    0.9870824  0.8417658
```

Apart from the Horseshoe prior, the ridge prior is also provided.

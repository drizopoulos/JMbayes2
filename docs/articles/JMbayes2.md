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
#> time: 12 sec
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
function from the **GLMMadaptive** package. In addition, the argument
`which_independent` can be used to determine which longitudinal outcomes
are to be assumed independent; here, as an illustration, we specify that
the first (i.e., serum bilirubin) and second (i.e., prothrombin time)
longitudinal outcomes are independent. To assume that all longitudinal
outcomes are independent, we can use
`jm(..., which_independent = "all")`. Because this joint model is more
complex, we increase the number of MCMC iterations, the number of
burn-in iterations, and the thinning per chain using the corresponding
control arguments:

``` r

fm2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)
fm3 <- mixed_model(ascites ~ year + sex, data = pbc2,
                   random = ~ 1 | id, family = binomial())

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
#> marginal    11742.42 12753.66 -6682.826
#> conditional 12810.75 12501.97 -6734.081
#> 
#> Random-effects covariance matrix:
#>                                           
#>        StdDev   Corr                      
#> (Intr) 0.9983 (Intr)   year (Intr)    year
#> year   0.1803 0.4389                      
#> (Intr) 0.7670                             
#> year   0.3194               -0.0228       
#> (Intr) 2.6885 0.5421 0.7130 0.2593  0.1190
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale            -0.3976 0.3138 -0.9951 0.2177 0.2170 1.0031
#> value(log(serBilir))  0.8589 0.1357  0.5782 1.1216 0.0000 1.0261
#> value(prothrombin)    0.1240 0.0559  0.0081 0.2261 0.0377 1.0047
#> value(ascites)        0.3290 0.0949  0.1607 0.5326 0.0000 1.0771
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)     0.7295 0.1703  0.3903  1.0594 0.0000 0.9997
#> year            0.2735 0.0368  0.2029  0.3464 0.0000 1.0006
#> sexfemale      -0.2701 0.1809 -0.6241  0.0899 0.1323 0.9996
#> year:sexfemale -0.0881 0.0391 -0.1667 -0.0135 0.0193 1.0017
#> sigma           0.3488 0.0067  0.3356  0.3623 0.0000 1.0026
#> 
#> Longitudinal outcome: prothrombin (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    10.9791 0.1713 10.6401 11.3211 0.0000 1.0005
#> year            0.2234 0.0764  0.0755  0.3766 0.0030 0.9998
#> sexfemale      -0.4324 0.1824 -0.7938 -0.0720 0.0177 1.0001
#> year:sexfemale  0.0290 0.0805 -0.1292  0.1859 0.7107 1.0001
#> sigma           1.0591 0.0205  1.0195  1.1002 0.0000 1.0034
#> 
#> Longitudinal outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept) -4.2309 0.6197 -5.4615 -3.0728 0.0000 1.0065
#> year         0.4519 0.0530  0.3491  0.5548 0.0000 1.0106
#> sexfemale   -0.6088 0.6034 -1.8064  0.5323 0.3113 0.9997
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 12000 
#> burn-in per chain: 2000 
#> thinning: 5 
#> time: 1.1 min
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
  "log(serBilir)" = ~ slope(log(serBilir)):sex,
  "prothrombin" = ~ area(prothrombin)
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
#> marginal    11787.72 12918.39 -6649.515
#> conditional 12723.89 12419.93 -6662.220
#> 
#> Random-effects covariance matrix:
#>                                           
#>        StdDev   Corr                      
#> (Intr) 0.9798 (Intr)   year (Intr)    year
#> year   0.1877 0.5560                      
#> (Intr) 0.7679                             
#> year   0.3203               -0.0308       
#> (Intr) 2.7220 0.6462 0.7503 0.2068  0.1126
#> 
#> Survival outcome:
#>                                   Mean  StDev    2.5%   97.5%      P   Rhat
#> sexfemale                      -0.0533 0.8458 -1.5114  1.8152 0.8917 1.0094
#> slope(log(serBilir)):sexmale    7.4865 2.4955  3.1349 13.1678 0.0010 1.0300
#> slope(log(serBilir)):sexfemale  6.0769 1.8765  2.1844  9.8518 0.0050 1.0204
#> area(prothrombin)               0.3136 0.1234  0.0652  0.5466 0.0160 1.0022
#> value(ascites)                  0.4818 0.1411  0.2035  0.7653 0.0003 1.0137
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)     0.7197 0.1677  0.3916  1.0445 0.0000 1.0006
#> year            0.2804 0.0369  0.2089  0.3542 0.0000 1.0016
#> sexfemale      -0.2669 0.1775 -0.6153  0.0811 0.1380 1.0006
#> year:sexfemale -0.0795 0.0388 -0.1564 -0.0030 0.0417 1.0006
#> sigma           0.3515 0.0069  0.3381  0.3655 0.0000 1.0088
#> 
#> Longitudinal outcome: prothrombin (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    10.9725 0.1707 10.6345 11.3092 0.0000 1.0001
#> year            0.2263 0.0764  0.0718  0.3759 0.0057 1.0029
#> sexfemale      -0.4237 0.1801 -0.7804 -0.0683 0.0197 0.9999
#> year:sexfemale  0.0257 0.0803 -0.1333  0.1891 0.7490 1.0021
#> sigma           1.0585 0.0206  1.0183  1.0996 0.0000 1.0023
#> 
#> Longitudinal outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept) -4.1845 0.6097 -5.4423 -3.0372 0.0000 1.0034
#> year         0.4750 0.0522  0.3784  0.5798 0.0000 1.0135
#> sexfemale   -0.6252 0.6029 -1.8127  0.5634 0.2903 1.0013
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
#> gammas.Mean                                 -0.05326779 -0.1707259
#> association.slope(log(serBilir)):sexmale     7.48652502  6.2842336
#> association.slope(log(serBilir)):sexfemale   6.07691289  5.2163420
#> association.area(prothrombin)                0.31359397  0.2350919
#> association.value(ascites)                   0.48177397  0.5454002
```

Apart from the Horseshoe prior, the ridge prior is also provided.

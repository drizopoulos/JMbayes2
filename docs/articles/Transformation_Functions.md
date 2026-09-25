# Transformation Functions for Functional Forms

## Functional Forms

### Simplified syntax

We have [previously
seen](https://drizopoulos.github.io/JMbayes2/articles/JMbayes2.html#functional-forms-1)
that function
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md) via its
`functional_forms` argument allows the specification of different
functional forms to link the longitudinal and event time outcomes. This
argument accepts either a single formula or a list of formulas per
longitudinal outcome with the terms we wish to include.

We will illustrate some of these possibilities using the PBC dataset. We
start by fitting a Cox model for the composite event transplantation or
death, including sex as a baseline covariate:

``` r

pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
```

Our aim is to assess the strength of the association between the risk of
the composite event and whether the patients experienced hepatomegaly
during follow-up. We will describe the patient-specific profiles over
time for this biomarker using a mixed-effects logistic model, where we
include an intercept and the time effect in both fixed and random
effects. The syntax to fit this model with
[`mixed_model()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
is:

``` r

fm <- mixed_model(hepatomegaly ~ year, data = pbc2, random = ~ year | id, 
                  family = binomial())
```

The default call to
[`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md) adds
the subject-specific linear predictor of the mixed-effects logistic
regression as a time-varying covariate in the survival relative risk
model:

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
#>   hepatomegaly: 1884
#> 
#>                  DIC     WAIC      LPML
#> marginal    3271.158 3264.337 -1632.730
#> conditional 4783.809 4846.361 -2664.249
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 3.3077 (Intr) 
#> year   0.5354 -0.1638
#> 
#> Survival outcome:
#>                        Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale           -0.4619 0.3028 -1.0504 0.1607 0.1287 1.0023
#> value(hepatomegaly)  0.3560 0.0476  0.2625 0.4475 0.0000 1.0394
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>               Mean  StDev    2.5%  97.5%     P   Rhat
#> (Intercept) 0.0514 0.2266 -0.3827 0.5064 0.826 1.0017
#> year        0.2871 0.0652  0.1640 0.4164 0.000 1.0017
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 15 sec
```

In the output, this is named `value(hepatomegaly)` to denote that the
current value functional form is used. That is, we assume that the risk
at a specific time t is associated with the value of the linear
predictor of the longitudinal outcome at the same time point t. In this
case, the subject-specific linear predictor denotes the log odds of
experiencing hepatomegaly at time t.

### Transformation functions

The fact that the default version of the current value functional form
is on the linear predictor scale of the mixed model may be problematic
to interpret when this linear predictor is connected with a nonlinear
link function to the mean of the longitudinal outcome. In these
situations, we may want to transform the subject-specific linear
predictor back to the scale of the outcome. To achieve this, we can use
a transformation function. Continuing on the previous example, we update
`jointFit1` by now linking the expit transformation of the linear
predictor (i.e., {\sf expit}(x) = \exp(x) / \\1 + \exp(x)\\) with the
risk of an event. This is done using the
[`vexpit()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md)
function:

``` r

jointFit2 <- update(jointFit1, functional_forms = ~ vexpit(value(hepatomegaly)))
summary(jointFit2)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = CoxFit, Mixed_objects = fm, time_var = "year", 
#>     functional_forms = ~vexpit(value(hepatomegaly)))
#> 
#> Data Descriptives:
#> Number of groups: 312        Number of events: 169 (54.2%)
#> Number of observations:
#>   hepatomegaly: 1884
#> 
#>                  DIC         WAIC          LPML
#> marginal    10948427 2.096815e+12 -211893560.62
#> conditional -2398343 4.898106e+03      -2681.99
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 3.3446 (Intr) 
#> year   0.5290 -0.3215
#> 
#> Survival outcome:
#>                                Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale                   -0.3525 0.2759 -0.9036 0.1922 0.2056 1.0043
#> vexpit(value(hepatomegaly))  3.2610 0.4417  2.4286 4.1611 0.0000 1.0220
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>               Mean  StDev    2.5%  97.5%      P   Rhat
#> (Intercept) 0.0736 0.2312 -0.3696 0.5419 0.7462 1.0046
#> year        0.2364 0.0627  0.1159 0.3595 0.0000 1.0059
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 16 sec
```

Other available functions to use in the definition of the
`functional_forms` argument are
[`vexp()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md) to
calculate the exponent,
[`vlog()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md) to
calculate the natural logarithm, and
[`vsqrt()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md) to
calculate the square root.

If we want to include the time-varying slope of the transformed linear
predictor, we also have the
[`Dexpit()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md) and
[`Dexp()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md)
functions available. As an example, we extend `jointFit2` by including
the derivative of the {\sf expit}() transformation:

``` r

forms <- ~ vexpit(value(hepatomegaly)) + Dexpit(slope(hepatomegaly))
jointFit3 <- update(jointFit1, functional_forms = forms)
summary(jointFit3)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = CoxFit, Mixed_objects = fm, time_var = "year", 
#>     functional_forms = forms)
#> 
#> Data Descriptives:
#> Number of groups: 312        Number of events: 169 (54.2%)
#> Number of observations:
#>   hepatomegaly: 1884
#> 
#>                       DIC         WAIC          LPML
#> marginal     3.710210e+22 3.101078e+47 -6.213371e+25
#> conditional -1.123152e+10 4.879972e+03 -2.707586e+03
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 3.4002 (Intr) 
#> year   0.5488 -0.4588
#> 
#> Survival outcome:
#>                                                    Mean  StDev    2.5%   97.5%
#> sexfemale                                       -0.3313 0.2846 -0.8546  0.2500
#> vexpit(value(hepatomegaly))                      3.4213 0.5315  2.4733  4.5506
#> Dexpit(value(hepatomegaly)):slope(hepatomegaly) -0.9514 0.4716 -1.9941 -0.1333
#>                                                      P   Rhat
#> sexfemale                                       0.2467 1.0025
#> vexpit(value(hepatomegaly))                     0.0000 1.0605
#> Dexpit(value(hepatomegaly)):slope(hepatomegaly) 0.0222 1.0332
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>               Mean  StDev    2.5%  97.5%      P   Rhat
#> (Intercept) 0.1099 0.2356 -0.3458 0.5814 0.6489 1.0062
#> year        0.1766 0.0639  0.0541 0.3055 0.0040 1.0278
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 18 sec
```

The call to `Dexpit(slope(hepatomegaly))` is internally transformed to
`Dexpit(value(hepatomegaly)):slope(hepatomegaly)`, which calculates the
derivative of the {\sf expit}() evaluated at the linear predictor times
the derivative of the linear predictor. This is because \frac{d}{dt}
{\sf expit}\\\eta(t)\\ = {\sf expit}\\\eta(t)\\ \\ \Bigl \[ 1 - {\sf
expit}\\\eta(t)\\ \Bigr \] \times \frac{d}{dt}\eta(t)

### The Slope functional form

As we have seen in previous examples, the
[`slope()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md)
function is used to specify the slope functional form d \eta(t)/dt.
According to the [definition of the
derivative](https://en.wikipedia.org/wiki/Derivative#Rigorous_definition),
this corresponds to the change in the longitudinal profile \\\eta(t +
\varepsilon) - \eta(t)\\/ \varepsilon as \varepsilon approaches zero.
However, the interpretation of this term may be challenging in some
settings. A possible alternative would be to increase the value of
\varepsilon, e.g., \varepsilon = 1. For example, if the time scale is
years, this would quantify the change of the longitudinal profile in the
last year before t.

To fit a joint model with such a term, we can use the `eps` and
`direction` arguments of the
[`slope()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md)
function. We illustrate this in the following example, in which we use
the serum bilirubin

``` r

gm <- lme(log(serBilir) ~ ns(year, 2), data = pbc2, random = ~ ns(year, 2) | id,
          control = lmeControl(opt = "optim"))
```

We first fit the joint model with time-varying slope term:

``` r

jFit1 <- jm(CoxFit, gm, time_var = "year",
            functional_forms = ~ value(log(serBilir)) + slope(log(serBilir)))
```

To specify that we want to include the change in the log serum bilirubin
levels during the last year before t, we specify `eps = 1` and
`direction = "back"` in the call to
[`slope()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md).
This calculates the term \\\eta(t) - \eta(t - \varepsilon)\\ /
\varepsilon for \varepsilon set equal to `eps = 1`:

``` r

jFit2 <- jm(CoxFit, gm, time_var = "year",
            functional_forms = ~ value(log(serBilir)) + 
              slope(log(serBilir), eps = 1, direction = "back"))
```

We compare the two fits

``` r

summary(jFit1)$Survival
#>                           Mean     StDev       2.5%     97.5%         P
#> sexfemale            -0.175050 0.2726780 -0.6789224 0.3783756 0.5268889
#> value(log(serBilir))  1.220027 0.1090261  1.0128844 1.4400938 0.0000000
#> slope(log(serBilir))  2.944295 0.6587145  1.7427147 4.2971061 0.0000000
#>                          Rhat
#> sexfemale            1.001977
#> value(log(serBilir)) 1.024068
#> slope(log(serBilir)) 1.027657

summary(jFit2)$Survival
#>                                                         Mean     StDev
#> sexfemale                                         -0.1774053 0.2716258
#> value(log(serBilir))                               1.2071915 0.1036586
#> slope(log(serBilir), eps = 1, direction = "back")  2.8198207 0.6105925
#>                                                         2.5%     97.5%
#> sexfemale                                         -0.6966948 0.3748471
#> value(log(serBilir))                               1.0068875 1.4087023
#> slope(log(serBilir), eps = 1, direction = "back")  1.6957148 4.1260694
#>                                                           P     Rhat
#> sexfemale                                         0.5366667 1.002943
#> value(log(serBilir))                              0.0000000 1.049348
#> slope(log(serBilir), eps = 1, direction = "back") 0.0000000 1.010129
```

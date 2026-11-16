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
#> marginal    3078.721 3068.122 -1535.149
#> conditional 4991.883 4851.605 -2659.798
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 3.2972 (Intr) 
#> year   0.5369 -0.1629
#> 
#> Survival outcome:
#>                        Mean StDev    2.5%  97.5%      P   Rhat
#> sexfemale           -0.4714 0.299 -1.0511 0.1216 0.1171 1.0012
#> value(hepatomegaly)  0.3521 0.049  0.2629 0.4554 0.0000 1.0265
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>               Mean  StDev    2.5%  97.5%      P   Rhat
#> (Intercept) 0.0527 0.2262 -0.3858 0.5010 0.8253 1.0043
#> year        0.2853 0.0615  0.1670 0.4072 0.0000 1.0219
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
#>                  DIC     WAIC      LPML
#> marginal    3065.389 3058.731 -1530.035
#> conditional 5022.700 4924.262 -2726.802
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 3.4064 (Intr) 
#> year   0.5426 -0.3475
#> 
#> Survival outcome:
#>                                Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale                   -0.3586 0.2816 -0.8643 0.2012 0.2031 1.0005
#> vexpit(value(hepatomegaly))  3.2322 0.4437  2.4012 4.1573 0.0000 1.0166
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>               Mean  StDev    2.5%  97.5%      P   Rhat
#> (Intercept) 0.0698 0.2347 -0.3758 0.5339 0.7647 1.0034
#> year        0.2359 0.0631  0.1193 0.3644 0.0000 1.0164
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
#>                  DIC     WAIC      LPML
#> marginal    3066.994 3058.858 -1530.209
#> conditional 4997.713 4894.311 -2713.970
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 3.4536 (Intr) 
#> year   0.5612 -0.4739
#> 
#> Survival outcome:
#>                                                    Mean  StDev    2.5%   97.5%
#> sexfemale                                       -0.3289 0.2898 -0.8727  0.2508
#> vexpit(value(hepatomegaly))                      3.4277 0.5446  2.4731  4.6946
#> Dexpit(value(hepatomegaly)):slope(hepatomegaly) -0.9622 0.4801 -2.1282 -0.1991
#>                                                      P   Rhat
#> sexfemale                                       0.2678 1.0040
#> vexpit(value(hepatomegaly))                     0.0000 1.2764
#> Dexpit(value(hepatomegaly)):slope(hepatomegaly) 0.0222 1.0466
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>               Mean  StDev    2.5%  97.5%      P   Rhat
#> (Intercept) 0.1206 0.2397 -0.3446 0.6016 0.6164 1.0178
#> year        0.1725 0.0665  0.0469 0.3056 0.0087 1.0548
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
#>                            Mean     StDev      2.5%    97.5%         P     Rhat
#> sexfemale            -0.1819164 0.2793902 -0.695733 0.395819 0.5155556 1.004030
#> value(log(serBilir))  1.2282014 0.1088604  1.009231 1.440659 0.0000000 1.092506
#> slope(log(serBilir))  2.9379169 0.6423208  1.732051 4.214884 0.0000000 1.032774

summary(jFit2)$Survival
#>                                                         Mean    StDev
#> sexfemale                                         -0.1738989 0.278299
#> value(log(serBilir))                               1.2083299 0.107649
#> slope(log(serBilir), eps = 1, direction = "back")  2.8098669 0.612553
#>                                                         2.5%    97.5%         P
#> sexfemale                                         -0.7011899 0.402321 0.5282222
#> value(log(serBilir))                               0.9996800 1.423360 0.0000000
#> slope(log(serBilir), eps = 1, direction = "back")  1.6327652 4.072773 0.0000000
#>                                                       Rhat
#> sexfemale                                         1.000880
#> value(log(serBilir))                              1.076197
#> slope(log(serBilir), eps = 1, direction = "back") 1.032615
```

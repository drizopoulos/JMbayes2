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
#> Number of Groups: 312        Number of events: 169 (54.2%)
#> Number of Observations:
#>   hepatomegaly: 1884
#> 
#>                  DIC     WAIC      LPML
#> marginal    3077.880 3066.424 -1533.651
#> conditional 4988.864 4850.268 -2651.221
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 3.2996 (Intr) 
#> year   0.5368 -0.1730
#> 
#> Survival Outcome:
#>                        Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale           -0.4733 0.3019 -1.0613 0.1357 0.1191 1.0011
#> value(hepatomegaly)  0.3591 0.0531  0.2641 0.4690 0.0000 1.0712
#> 
#> Longitudinal Outcome: hepatomegaly (family = binomial, link = logit)
#>               Mean  StDev    2.5%  97.5%      P   Rhat
#> (Intercept) 0.0527 0.2282 -0.3928 0.5037 0.8244 1.0037
#> year        0.2836 0.0646  0.1561 0.4083 0.0000 1.0248
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 16 sec
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
predictor (i.e., \mbox{expit}(x) = \exp(x) / \\1 + \exp(x)\\) with the
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
#> Number of Groups: 312        Number of events: 169 (54.2%)
#> Number of Observations:
#>   hepatomegaly: 1884
#> 
#>                  DIC     WAIC      LPML
#> marginal    3068.775 3059.436 -1530.111
#> conditional 5005.254 4904.042 -2701.585
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 3.3268 (Intr) 
#> year   0.5302 -0.3244
#> 
#> Survival Outcome:
#>                                Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale                   -0.3479 0.2838 -0.8972 0.2201 0.2196 1.0019
#> vexpit(value(hepatomegaly))  3.2874 0.4691  2.3812 4.2400 0.0000 1.0166
#> 
#> Longitudinal Outcome: hepatomegaly (family = binomial, link = logit)
#>               Mean  StDev    2.5%  97.5%      P   Rhat
#> (Intercept) 0.0657 0.2354 -0.3840 0.5423 0.7862 1.0015
#> year        0.2377 0.0619  0.1191 0.3605 0.0004 1.0017
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 18 sec
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
the derivative of the \mbox{expit}() transformation:

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
#> Number of Groups: 312        Number of events: 169 (54.2%)
#> Number of Observations:
#>   hepatomegaly: 1884
#> 
#>                  DIC     WAIC      LPML
#> marginal    3067.439 3058.267 -1529.488
#> conditional 4990.422 4886.822 -2701.964
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 3.4301 (Intr) 
#> year   0.5544 -0.4777
#> 
#> Survival Outcome:
#>                                                    Mean  StDev    2.5%   97.5%
#> sexfemale                                       -0.3234 0.2953 -0.8756  0.2695
#> vexpit(value(hepatomegaly))                      3.2922 0.5259  2.3173  4.3695
#> Dexpit(value(hepatomegaly)):slope(hepatomegaly) -0.9686 0.4667 -2.0680 -0.2517
#>                                                      P   Rhat
#> sexfemale                                       0.2764 1.0079
#> vexpit(value(hepatomegaly))                     0.0000 1.0142
#> Dexpit(value(hepatomegaly)):slope(hepatomegaly) 0.0071 1.3755
#> 
#> Longitudinal Outcome: hepatomegaly (family = binomial, link = logit)
#>               Mean  StDev    2.5%  97.5%      P   Rhat
#> (Intercept) 0.1410 0.2439 -0.3245 0.6345 0.5656 1.0179
#> year        0.1603 0.0688  0.0262 0.2970 0.0156 1.1988
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 22 sec
```

The call to `Dexpit(slope(hepatomegaly))` is internally transformed to
`Dexpit(value(hepatomegaly)):slope(hepatomegaly)`, which calculates the
derivative of the \mbox{expit}() evaluated at the linear predictor times
the derivative of the linear predictor. This is because \frac{d}{dt}
\mbox{expit}\\\eta(t)\\ = \mbox{expit}\\\eta(t)\\ \\ \Bigl \[ 1 -
\mbox{expit}\\\eta(t)\\ \Bigr \] \times \frac{d}{dt}\eta(t)

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
#>                            Mean     StDev       2.5%    97.5%         P
#> sexfemale            -0.1683743 0.2795170 -0.6744899 0.405122 0.5535556
#> value(log(serBilir))  1.2223588 0.1054294  1.0243374 1.429878 0.0000000
#> slope(log(serBilir))  2.9624164 0.6749065  1.7579233 4.393638 0.0000000
#>                          Rhat
#> sexfemale            1.008869
#> value(log(serBilir)) 1.010025
#> slope(log(serBilir)) 1.043587

summary(jFit2)$Survival
#>                                                         Mean     StDev
#> sexfemale                                         -0.1647358 0.2732504
#> value(log(serBilir))                               1.2098188 0.1075358
#> slope(log(serBilir), eps = 1, direction = "back")  2.8433863 0.6355504
#>                                                         2.5%     97.5%
#> sexfemale                                         -0.6761072 0.3902834
#> value(log(serBilir))                               0.9971623 1.4245781
#> slope(log(serBilir), eps = 1, direction = "back")  1.7099160 4.1559186
#>                                                           P     Rhat
#> sexfemale                                         0.5344444 1.007342
#> value(log(serBilir))                              0.0000000 1.013989
#> slope(log(serBilir), eps = 1, direction = "back") 0.0000000 1.044426
```

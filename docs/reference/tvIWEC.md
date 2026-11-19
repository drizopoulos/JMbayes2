# Integrated Weighted Expected Cost for Dynamic Predictions under Left-Truncated and Right-Censored Data

Computes estimates of the integrated weighted expected cost (IWEC) for a
given weight function within the beta family at a specified horizon time
(`Thoriz`), using longitudinal information available up to a given
starting time point (`Tstart`). The function supports joint models, Cox
proportional hazards models and predictions obtained from any method. It
provides inverse probability weighting and model-based estimation
methods.

## Usage

``` r
tvIWEC(object, ...)

# Default S3 method
tvIWEC(object, Time, event, LT_Time = NULL, Tstart,
       shape1 = NULL, shape2 = NULL, mod = NULL, interv = NULL,
       mass = 0.95, Thoriz = NULL, Dt = NULL,
       type_weights = c("model-based", "IPCW"), preds_MB = NULL,
       opt_interv = c(1, 500), name_model = NULL, ...)

# S3 method for class 'jm'
tvIWEC(object, newdata, Tstart,
       shape1 = NULL, shape2 = NULL, mod = NULL, interv = NULL, mass = 0.95,
       Thoriz = NULL, Dt = NULL,
       type_weights = c("model-based", "IPCW"),
       opt_interv = c(1, 500), name_model = NULL, ...)

# S3 method for class 'coxph'
tvIWEC(object, newdata, Tstart,
       shape1 = NULL, shape2 = NULL, mod = NULL, interv = NULL, mass = 0.95,
       Thoriz = NULL, Dt = NULL,
       type_weights = c("model-based", "IPCW"),
       opt_interv = c(1, 500), name_model = NULL, ...)
```

## Arguments

- object:

  an object inheriting from class `jm` or `coxph`, or a numeric vector
  containing the dynamic predictions to be evaluated.

- newdata:

  a data.frame that contains the longitudinal and covariate information
  for the subjects for which prediction of survival probabilities is
  required. If a `jm` object is used, the names of the variables in this
  data.frame must be the same as in the data.frames that were used to
  fit the linear mixed effects and the event process model that were
  supplied as the two first argument of
  [`jm`](https://drizopoulos.github.io/JMbayes2/reference/jm.md). If a
  `coxph` is used it only applies to the event process. If the default
  method is used, the argument is not needed.

- Tstart:

  numeric scalar denoting the time point up to which longitudinal
  information is to be used to derive predictions.

- Time, event, LT_Time:

  numeric vectors that must be specified when using the default method
  to compute the metric by directly using a vector of predictions,
  observed times (`Time`), the event indicator (`event`) and (possibly)
  left-truncated times (`LT_Time`).

- shape1, shape2:

  numeric scalars greater or equal to 0 denoting the shape parameteres
  for the beta family of weight functions (see below for details).

- mod, mass, interv:

  two numeric scalars between 0 and 1 and a numeric vector with two
  scalars between 0 and 1, denoting the mode of the weight function and
  the mass (default is 95%) to be allocated within the interval defined
  by the two values in the numeric vector. These arguments can be used
  to compute the shape parameters. See below for details.

- Thoriz:

  numeric scalar denoting the time point for which a prediction of the
  survival status is of interest; `Thoriz` must be later than `Tstart`
  and either `Dt` or `Thoriz` must be specified. If `Thoriz` is `NULL`
  is set equal to `Tstart + Dt`.

- Dt:

  numeric scalar denoting the length of the time interval of prediction;
  either `Dt` or `Thoriz` must be specified.

- type_weights:

  character string denoting the type of weights to use to account for
  censorting. Options are model-based (default) and inverse probability
  weighting (using the Kaplan-Meier estimate of the censoring
  distribution, and adjusting for left truncation if necessary).

- preds_MB:

  a numeric vector containing the model-based weights to be used for the
  individuals censored within (`Tstart`, `Thoriz`\]. It must be
  specified when using the default method and the `model-based`
  estimator.

- opt_interv:

  a numeric vector containing the end-points of the interval to be
  searched for the root. This argumentis used when `mod`, `mass` and
  `interv` are provided, and shape parameters are computed (see details
  below).

- name_model:

  character string denoting the name of the model or object used for
  prediction evaluation; default is `NULL`, meaning that the name of the
  object in the environment will be used.

- ...:

  additional arguments passed to
  [`predict.jm()`](https://drizopoulos.github.io/JMbayes2/reference/predict.md).

## Details

The IWEC is implemented using the beta family of weight functions, which
is defined by the two-parameter form \\w(c;\text{shape1},
\text{shape2})=c^{\text{shape1}-1}(1-c)^{\text{shape2}-1}/B(\text{shape1},
\text{shape2})\\, where \\\text{shape1}, \text{shape2}\geq 0\\ are shape
parameters and \\B(\text{shape1},\text{shape2})=1\\ if
\\\text{shape1}=\text{shape2}=0\\, and otherwise denotes the beta
function.

When `mod`, `mass` and `interv` are provided we compute the shape
parameters `shape1` and `shape2` by solving the system consisting of the
two equations: \\\text{mod}=\frac{\text{shape1} - 1}{\text{shape1} +
\text{shape2} - 2}\\ and \\\text{Pr}(X \in \text{interv}) =
\text{mass}\\, where \\X\\ follows a beta distribution with shape
parameters `shape1` and `shape2`. To solve this system, the
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) function is
internally used. The user can change the `interval` argument of the
[`uniroot()`](https://rdrr.io/r/stats/uniroot.html) function using the
`opt_interv` argument.

## Value

A list of class `tvMetric` with components:

- metric:

  a numeric scalar denoting the estimated IWEC.

- metric_type:

  the character string `'IWEC'`.

- nr:

  a numeric scalar denoting the number of subjects at risk at time
  `Tstart`.

- nint:

  a numeric scalar denoting the number of subjects with an event in
  (`Tstart`, `Thoriz`\].

- ncens:

  a numeric scalar denoting the number of subjects with a censored time
  in (`Tstart`, `Thoriz`\].

- Tstart:

  a copy of the `Tstart` argument.

- Thoriz:

  a copy of the `Thoriz` argument.

- type_weights:

  a copy of the `type_weights` argument.

- metric_info:

  a copy of `shape1` and `shape2`.

- LT_MS:

  a logical that is `TRUE` if left-truncated times have been used.

- classObject:

  the class of `object`.

- nameObject:

  the name of the model; either a copy of the `name_model` argument, or
  the name of `object`.

## Author

Arnau Garcia-Fernandez <a.garciafernandez@erasmusmc.nl>, Pedro
Miranda-Afonso <p.mirandaafonso@erasmusmc.nl> and Dimitris Rizopoulos
<d.rizopoulos@erasmusmc.nl>

## See also

[`decision_curves`](https://drizopoulos.github.io/JMbayes2/reference/decision_curves.md),
[`tvNB`](https://drizopoulos.github.io/JMbayes2/reference/tvNB.md),
[`tvEC`](https://drizopoulos.github.io/JMbayes2/reference/tvEC.md),[`predict`](https://rdrr.io/r/stats/predict.html),
[`jm`](https://drizopoulos.github.io/JMbayes2/reference/jm.md)

## Examples

``` r
# \donttest{
# We fit a multivariate joint model
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
fm1 <- lme(log(serBilir) ~ ns(year, 3) * sex, data = pbc2,
           random = ~ ns(year, 3) | id, control = lmeControl(opt = 'optim'))
fm2 <- lme(prothrombin ~ ns(year, 2) * sex, data = pbc2,
           random = ~ ns(year, 2) | id, control = lmeControl(opt = 'optim'))

jointFit <- jm(CoxFit, list(fm1, fm2), time_var = "year", n_chains = 1L)

# IWEC with a weight function in which relevant decision thresholds are (0.1, 0.3),
# the mode is 0.2 and the mass within the interval is 80 percent.
# We use landmark time 6 and prediction horizon 4.
tvIWEC(jointFit, newdata = pbc2, Tstart = 6, Dt = 4, mod = 0.2,
       interv = c(0.1, 0.3), mass = 0.8, cores = 1L)
#> 
#> Accuracy Metric for the Joint Model 'jointFit'
#> 
#> Estimated Integrated Weighted Expected Cost
#>  with shape parameters 6.1 and 21.2: 0.0884
#> At time: 10
#> For the 166 subjects at risk at time 6
#> Number of subjects with an event in [6, 10): 33
#> Number of subjects with a censored time in [6, 10): 82
#> Accounting for censoring using model-based weights
#> 

# Considering the Cox model, the IPW approach for estimation and
# shape parameters 12 and 38
tvIWEC(CoxFit, newdata = pbc2.id, Tstart = 6, Dt = 4, shape1 = 12, shape2 = 38,
       type_weights = "IPCW")
#> 
#> Accuracy Metric for the Cox Model 'CoxFit'
#> 
#> Estimated Integrated Weighted Expected Cost
#>  with shape parameters 12.0 and 38.0: 0.1596
#> At time: 10
#> For the 166 subjects at risk at time 6
#> Number of subjects with an event in [6, 10): 42
#> Number of subjects with a censored time in [6, 10): 73
#> Accounting for censoring using 
#> inverse probability of censoring Kaplan-Meier weights
#> 

# Considering the vectors with predictions, observed times, event indicator and
# left-truncated times. We use the default method:
pbc2$lt <- unlist(by(pbc2, pbc2$id, function(x) rep(runif(1, 0, x$years), nrow(x))))
pbc2.id$lt <- runif(nrow(pbc2.id), 0, pbc2.id$years)
sfit0 <- summary(survfit(CoxFit, newdata = pbc2.id), times = 5)
sfit <- summary(survfit(CoxFit, newdata = pbc2.id), times = 5 + 3)
predictions <- 1.0 - (as.matrix(sfit$surv)[1L, ] / as.matrix(sfit0$surv)[1L, ])
tvIWEC(predictions, pbc2.id$years, pbc2.id$status2, pbc2.id$lt, Tstart = 5, Dt = 3,
       type_weights = "IPCW", shape1 = 12, shape2 = 35)
#> 
#> Accuracy Metric for the predictions 'predictions',
#>  observed times 'pbc2.id$years', event indicator 'pbc2.id$status2'
#>  and left-truncated times 'pbc2.id$lt'
#> 
#> Estimated Integrated Weighted Expected Cost
#>  with shape parameters 12.0 and 35.0: 0.1748
#> At time: 8
#> For the 131 subjects at risk at time 5
#> Number of subjects with an event in [5, 8): 31
#> Number of subjects with a censored time in [5, 8): 45
#> Accounting for left truncation and right censoring using 
#>  nonparametric inverse probability weights
#> 

# Now, we consider the default method, but using the model-based estimator:
ind_censoring <- pbc2.id$years > 5 & pbc2.id$lt < 5 & pbc2.id$years < 8 & pbc2.id$status2 == 0
sfit2_8 <- summary(survfit(CoxFit, newdata = pbc2.id[ind_censoring, ]),
                 times = 8)
sfit_cens <- summary(survfit(CoxFit, newdata = pbc2.id[ind_censoring, ]),
                     times = pbc2.id$years[ind_censoring])
preds_MB <- rep(NA, length(pbc2.id$years))
preds_MB[ind_censoring] <- 1.0 - (as.matrix(sfit2_8$surv)[1L, ] / as.matrix(sfit_cens$surv)[1L, ])
tvIWEC(predictions, Time = pbc2.id$years, event = pbc2.id$status2, LT_Time = pbc2.id$lt, 
       Tstart = 5, Dt = 3, type_weights = "model-based", preds_MB = preds_MB, 
       mod = 0.33, interv = c(0.125, 0.46), mass = 0.7)
#> 
#> Accuracy Metric for the predictions 'predictions',
#>  observed times 'pbc2.id$years', event indicator 'pbc2.id$status2'
#>  and left-truncated times 'pbc2.id$lt'
#> 
#> Estimated Integrated Weighted Expected Cost
#>  with shape parameters 3.5 and 6.0: 0.1368
#> At time: 8
#> For the 131 subjects at risk at time 5
#> Number of subjects with an event in [5, 8): 31
#> Number of subjects with a censored time in [5, 8): 45
#> Accounting for left truncation and right censoring using 
#>  model-based weights
#> 
# }
```

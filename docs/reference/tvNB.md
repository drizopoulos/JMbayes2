# Net Benefit for Dynamic Predictions under Left-Truncated and Right-Censored Data

Computes estimates of the net benefit for a given decision threshold
(`decision_threshold`) at a specified horizon time (`Thoriz`), using
longitudinal information available up to a given starting time point
(`Tstart`). The function supports joint models, Cox proportional hazards
models and predictions obtained from any method. It provides inverse
probability weighting and model-based estimation methods.

## Usage

``` r
tvNB(object, ...)

# Default S3 method
tvNB(object, Time, event, LT_Time = NULL, Tstart,
     decision_threshold = NULL, Thoriz = NULL, Dt = NULL,
     type_weights = c("model-based", "IPCW"),
     preds_MB = NULL, name_model = NULL, ...)

# S3 method for class 'jm'
tvNB(object, newdata, Tstart,
     decision_threshold = NULL, Thoriz = NULL, Dt = NULL,
     type_weights = c("model-based", "IPCW"), name_model = NULL, ...)

# S3 method for class 'coxph'
tvNB(object, newdata, Tstart,
     decision_threshold = NULL, Thoriz = NULL, Dt = NULL,
     type_weights = c("model-based", "IPCW"), name_model = NULL, ...)
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

- decision_threshold:

  numeric scalar between 0 and 1 denoting the decision threshold used in
  the net benefit.

- Time, event, LT_Time:

  numeric vectors that must be specified when using the default method
  to compute the metric by directly using a vector of predictions,
  observed times (`Time`), the event indicator (`event`) and (possibly)
  left-truncated times (`LT_Time`).

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

- name_model:

  character string denoting the name of the model or object used for
  prediction evaluation; default is `NULL`, meaning that the name of the
  object in the environment will be used.

- ...:

  additional arguments passed to
  [`predict.jm()`](https://drizopoulos.github.io/JMbayes2/reference/predict.md).

## Value

A list of class `tvMetric` with components:

- metric:

  a numeric scalar denoting the estimated net benefit.

- metric_type:

  the character string `'NB'`.

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

  a copy of the `decision_threshold` argument.

- LT_MS:

  a logical that is `TRUE` if left-truncated times have been used.

- classObject:

  the class of `object`.

- nameObject:

  the name of the model; either a copy of the `name_model` argument, or
  the name of `object`.

## References

Vickers, Andrew J and Elkin, Elena B (2006). Decision curve analysis: a
novel method for evaluating prediction models. *Medical Decision Making*
**26**, 565–574.

## Author

Arnau Garcia-Fernandez <a.garciafernandez@erasmusmc.nl>, Pedro
Miranda-Afonso <p.mirandaafonso@erasmusmc.nl> and Dimitris Rizopoulos
<d.rizopoulos@erasmusmc.nl>

## See also

[`decision_curves`](https://drizopoulos.github.io/JMbayes2/reference/decision_curves.md),
[`tvIWEC`](https://drizopoulos.github.io/JMbayes2/reference/tvIWEC.md),
[`tvEC`](https://drizopoulos.github.io/JMbayes2/reference/tvEC.md),
[`predict`](https://rdrr.io/r/stats/predict.html),
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

# Net benefit with decision threshold 0.3, prediction horizon 7
# and landmark time 4
tvNB(jointFit, newdata = pbc2, Tstart = 4, Dt = 3, decision_threshold = 0.3,
     cores = 1L)
#> 
#> Accuracy Metric for the Joint Model 'jointFit'
#> 
#> Estimated Net Benefit with decision threshold 0.3: 0.0892
#> At time: 7
#> For the 225 subjects at risk at time 4
#> Number of subjects with an event in [4, 7): 34
#> Number of subjects with a censored time in [4, 7): 62
#> Accounting for censoring using model-based weights
#> 

# Considering the Cox model and IPW approach for estimation:
tvNB(CoxFit, newdata = pbc2.id, Tstart = 4, Dt = 3, decision_threshold = 0.3,
     type_weights = "IPCW")
#> 
#> Accuracy Metric for the Cox Model 'CoxFit'
#> 
#> Estimated Net Benefit with decision threshold 0.3: -0.0016
#> At time: 7
#> For the 225 subjects at risk at time 4
#> Number of subjects with an event in [4, 7): 50
#> Number of subjects with a censored time in [4, 7): 46
#> Accounting for censoring using 
#> inverse probability of censoring Kaplan-Meier weights
#> 


# Considering the vectors with predictions, observed times, event indicator and
# left-truncated times. We use the default method:
pbc2.id$lt <- runif(nrow(pbc2.id), 0, pbc2.id$years)
sfit0 <- summary(survfit(CoxFit, newdata = pbc2.id), times = 5)
sfit <- summary(survfit(CoxFit, newdata = pbc2.id), times = 5 + 3)
predictions <- 1.0 - (as.matrix(sfit$surv)[1L, ] / as.matrix(sfit0$surv)[1L, ])
tvNB(predictions, pbc2.id$years, pbc2.id$status2, pbc2.id$lt, Tstart = 5, Thoriz = 8,
     type_weights = "IPCW", decision_threshold = 0.25)
#> 
#> Accuracy Metric for the predictions 'predictions',
#>  observed times 'pbc2.id$years', event indicator 'pbc2.id$status2'
#>  and left-truncated times 'pbc2.id$lt'
#> 
#> Estimated Net Benefit with decision threshold 0.25: 0.0261
#> At time: 8
#> For the 125 subjects at risk at time 5
#> Number of subjects with an event in [5, 8): 33
#> Number of subjects with a censored time in [5, 8): 44
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
tvNB(predictions, Time = pbc2.id$years, event = pbc2.id$status2, LT_Time = pbc2.id$lt, 
     Tstart = 5, Dt = 3, type_weights = "model-based", preds_MB = preds_MB, 
     decision_threshold = 0.3)
#> 
#> Accuracy Metric for the predictions 'predictions',
#>  observed times 'pbc2.id$years', event indicator 'pbc2.id$status2'
#>  and left-truncated times 'pbc2.id$lt'
#> 
#> Estimated Net Benefit with decision threshold 0.3: 0.0124
#> At time: 8
#> For the 125 subjects at risk at time 5
#> Number of subjects with an event in [5, 8): 33
#> Number of subjects with a censored time in [5, 8): 44
#> Accounting for left truncation and right censoring using 
#>  model-based weights
#> 
# }
```

# Time-Dependent Predictive Accuracy Measures for Joint Models

Using the available longitudinal information up to a starting time
point, these functions compute estimates of the ROC curve and the AUC,
the Brier score and expected predictive cross-entropy at a horizon time
point based on joint models.

## Usage

``` r
tvROC(object, newdata, Tstart, ...)

# S3 method for class 'jm'
tvROC(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, type_weights = c("model-based", "IPCW"), ...)

tvAUC(object, newdata, Tstart, ...)

# S3 method for class 'jm'
tvAUC(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, type_weights = c("model-based", "IPCW"), ...)

# S3 method for class 'tvROC'
tvAUC(object, ...)

calibration_plot(object, newdata, Tstart, ...)

# S3 method for class 'jm'
calibration_plot(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, df_ns = NULL, plot = TRUE,
    col = "red", lty = 1, lwd = 1,
    add_CI = TRUE, col_CI = "lightgrey",
    add_density = TRUE, col_dens = "grey",
    xlab = "Predicted Probabilities",
    ylab = "Observed Probabilities", main = "", ...)

calibration_metrics(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, df_ns = NULL, ...)

tvBrier(object, newdata, Tstart, ...)

# S3 method for class 'jm'
tvBrier(object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
    integrated = FALSE, type_weights = c("model-based", "IPCW"),
    model_weights = NULL, eventData_fun = NULL,
    parallel = c("snow", "multicore"),
    cores = parallelly::availableCores(omit = 1L), ...)

tvEPCE(object, newdata, Tstart, Thoriz = NULL, Dt = NULL, eps = 0.001,
    model_weights = NULL, eventData_fun = NULL,
    parallel = c("snow", "multicore"),
    cores = parallelly::availableCores(omit = 1L), ...)

create_folds(data, V = 5, id_var = "id",
    method = c("CV", "Bootstrap"), strata = NULL, seed = 123L)
```

## Arguments

- object:

  an object inheriting from class `jm`, except for `tvAUC.tvROC()` where
  this is an object of class `tvROC`. For `tvBrier()` and `tvEPCE()` it
  can also be a library of joint models.

- newdata:

  a data.frame that contains the longitudinal and covariate information
  for the subjects for which prediction of survival probabilities is
  required. The names of the variables in this data.frame must be the
  same as in the data.frames that were used to fit the linear mixed
  effects and the event process model that were supplied as the two
  first argument of
  [`jm`](https://drizopoulos.github.io/JMbayes2/reference/jm.md).

- Tstart:

  numeric scalar denoting the time point up to which longitudinal
  information is to be used to derive predictions.

- Thoriz:

  numeric scalar denoting the time point for which a prediction of the
  survival status is of interest; `Thoriz` must be later than `Tstart`
  and either `Dt` or `Thoriz` must be specified. If `Thoriz` is `NULL`
  is set equal to `Tstart + Dt`.

- Dt:

  numeric scalar denoting the length of the time interval of prediction;
  either `Dt` or `Thoriz` must be specified.

- integrated:

  logical; if `TRUE` the integrated Brier score is calculated.

- type_weights:

  character string denoting the type of weights to use to account for
  censorting. Options are model-based (default) and inverse probability
  of censoring weighting (using the Kaplan-Meier estimate of the
  censoring distribution).

- eps:

  numeric scalar used in the approximation of the hazard function.

- model_weights:

  a numeric vector of weights to combine predictions when `object` is a
  list of joint models of class `"jmList"`.

- eventData_fun:

  a function that takes as input the `newdata` and produces the dataset
  used for the event process model. This is useful when, for example,
  the event process model contains other time-varying covariates. It is
  important that this function does not alter the ordering of the
  subjects in `newdata`.

- parallel:

  character string; what type of parallel computing to use.

- cores:

  integer denoting the number of cores to be used when a library of
  joint models has been provided in `object`. If `cores = 1`, no
  parallel computing is used.

- df_ns:

  the degrees of freedom for the natural cubic spline of the cloglog
  transformation of the predicted probabilities used in the Cox model
  that assesses calibration. The default is 3 unless there are less than
  25 events in the interval (`Tstart`, `Thoriz`\] in which case it is 2.

- plot:

  logical; should a plot be produced. If `FALSE`, a list is returned
  with the observed and predicted probabilities.

- add_CI:

  logical; should 0.95 pointwise confidence intervals be added around
  the calibration line.

- col_CI:

  character; the color of the shaded area representing the 0.95
  pointwise confidence intervals around the calibration line.

- add_density:

  logical; should the kernal density estimation of the predicted
  probabilities be superimposed in the calibration plot.

- col, lwd, lty, col_dens, xlab, ylab, main:

  graphical parameters.

- data:

  the data.frame to split in folds.

- V:

  numeric scalar denoting the number of folds for cross-validation or
  the number of sample for the Bootstrap methods.

- id_var:

  character string denoting the name of the subject id variable in
  `data`.

- strata:

  character vector with the names of stratifying variables.

- method:

  character string indicating which method to use to create the training
  and testing datasets in `create_folds()`. The default is V-fold
  cross-validation. For the `Bootstrap` option, `V` samples with
  replacement from the original dataset are proruced as training data.
  The testing data contains the subjects that were not selected in the
  respective Bootstrap sample.

- seed:

  integer denoting the seed.

- ...:

  additional arguments passed to
  [`predict.jm()`](https://drizopoulos.github.io/JMbayes2/reference/predict.md).

## Value

A list of class `tvAUC` with components:

- auc:

  a numeric scalar denoting the estimated prediction error.

- Tstart:

  a copy of the `Tstart` argument.

- Thoriz:

  a copy of the `Thoriz` argument.

- nr:

  a numeric scalar denoting the number of subjects at risk at time
  `Tstart`.

- classObject:

  the class of `object`.

- nameObject:

  the name of `object`.

A list of class `tvROC` with components:

- TP, FP, nTP, nFN, nTN, qSN, qSP, qOverall:

  accuracy indexes.

- F1score, Youden:

  numeric scalars with the optimal cut-point using the F1 score and the
  Youden index.

- thr:

  numeric vector of thresholds.

- Tstart:

  a copy of the `Tstart` argument.

- Thoriz:

  a copy of the `Thoriz` argument.

- nr:

  a numeric scalar denoting the number of subjects at risk at time
  `Tstart`.

- classObject:

  the class of `object`.

- nameObject:

  the name of `object`.

## References

Antolini, L., Boracchi, P., and Biganzoli, E. (2005). A time-dependent
discrimination index for survival data. *Statistics in Medicine* **24**,
3927–3944.

Commenges, D., Liquet, B., and Proust-Lima, C. (2012). Choice of
prognostic estimators in joint models by estimating differences of
expected conditional Kullback-Leibler risks. *Biometrics* **68**,
380–387.

Harrell, F., Kerry, L. and Mark, D. (1996). Multivariable prognostic
models: issues in developing models, evaluating assumptions and
adequacy, and measuring and reducing errors. *Statistics in Medicine*
**15**, 361–387.

Heagerty, P. and Zheng, Y. (2005). Survival model predictive accuracy
and ROC curves. *Biometrics* **61**, 92–105.

Rizopoulos, D. (2016). The R package JMbayes for fitting joint models
for longitudinal and time-to-event data using MCMC. *Journal of
Statistical Software* **72(7)**, 1–45. doi:10.18637/jss.v072.i07.

Rizopoulos, D. (2012) *Joint Models for Longitudinal and Time-to-Event
Data: with Applications in R*. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in
joint models for longitudinal and time-to-event data. *Biometrics*
**67**, 819–829.

Rizopoulos, D., Molenberghs, G. and Lesaffre, E.M.E.H. (2017). Dynamic
predictions with time-dependent covariates in survival analysis using
joint modeling and landmarking. *Biometrical Journal* **59**, 1261–1276.

## Author

Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>

## See also

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
fm3 <- mixed_model(ascites ~ year * sex, data = pbc2,
                   random = ~ year | id, family = binomial())

jointFit <- jm(CoxFit, list(fm1, fm2, fm3), time_var = "year", n_chains = 1L)

roc <- tvROC(jointFit, newdata = pbc2, Tstart = 4, Dt = 3, cores = 1L)
roc
#> 
#>  Time-dependent Sensitivity and Specificity for the Joint Model jointFit
#> 
#> At time: 7
#> Using information up to time: 4 (225 subjects still at risk)
#> Accounting for censoring using model-based weights
#> 
#>    cut-off      SN     SP  
#> 1     0.00 0.00000 1.0000  
#> 2     0.01 0.01306 0.9978  
#> 3     0.03 0.02703 0.9958  
#> 4     0.05 0.04880 0.9958  
#> 5     0.06 0.06077 0.9933  
#> 6     0.07 0.06077 0.9877  
#> 7     0.08 0.08253 0.9877  
#> 8     0.09 0.10430 0.9877  
#> 9     0.10 0.12607 0.9877  
#> 10    0.13 0.14478 0.9813  
#> 11    0.17 0.15942 0.9795  
#> 12    0.18 0.18119 0.9795  
#> 13    0.19 0.18119 0.9739  
#> 14    0.20 0.23145 0.9700  
#> 15    0.21 0.28848 0.9679  
#> 16    0.22 0.30214 0.9658  
#> 17    0.23 0.32391 0.9602  
#> 18    0.25 0.34567 0.9602  
#> 19    0.26 0.36360 0.9481  
#> 20    0.30 0.36360 0.9425  
#> 21    0.32 0.38537 0.9425  
#> 22    0.34 0.40714 0.9425  
#> 23    0.36 0.42891 0.9369  
#> 24    0.37 0.45067 0.9369  
#> 25    0.41 0.49421 0.9369  
#> 26    0.45 0.50564 0.9287  
#> 27    0.46 0.50564 0.9231  
#> 28    0.47 0.51426 0.9141  
#> 29    0.51 0.52846 0.9066  
#> 30    0.54 0.53187 0.8963  
#> 31    0.55 0.55364 0.8851  
#> 32    0.56 0.56056 0.8813  
#> 33    0.57 0.60409 0.8758  
#> 34    0.58 0.64763 0.8758  
#> 35    0.60 0.66939 0.8758  
#> 36    0.61 0.66939 0.8702  
#> 37    0.62 0.66939 0.8590  
#> 38    0.63 0.69619 0.8547  
#> 39    0.65 0.71796 0.8547  
#> 40    0.66 0.71796 0.8491  
#> 41    0.68 0.71908 0.8382  
#> 42    0.69 0.74085 0.8382 *
#> 43    0.70 0.74085 0.8327  
#> 44    0.71 0.74530 0.8282  
#> 45    0.72 0.74530 0.8226  
#> 46    0.73 0.74530 0.8059  
#> 47    0.74 0.76707 0.8059  
#> 48    0.75 0.77172 0.7791  
#> 49    0.77 0.77398 0.7630  
#> 50    0.78 0.77398 0.7574  
#> 51    0.79 0.77398 0.7462  
#> 52    0.80 0.79806 0.7412  
#> 53    0.81 0.82369 0.7255  
#> 54    0.82 0.83206 0.6997  
#> 55    0.83 0.85792 0.6784  
#> 56    0.84 0.85843 0.6618  
#> 57    0.85 0.86268 0.6349  
#> 58    0.86 0.86854 0.5918  
#> 59    0.87 0.87327 0.5762  
#> 60    0.88 0.87327 0.5483  
#> 61    0.89 0.87385 0.5373  
#> 62    0.90 0.88081 0.4609  
#> 63    0.91 0.90438 0.4390  
#> 64    0.92 0.92843 0.4228  
#> 65    0.93 0.92945 0.3784  
#> 66    0.94 0.93024 0.3339  
#> 67    0.95 0.99618 0.2783  
#> 68    0.96 0.99936 0.1897  
#> 69    0.97 0.99963 0.1172  
#> 70    0.98 0.99990 0.0279  
#> 71    0.99 1.00000 0.0000  
#> 
tvAUC(roc)
#> 
#>  Time-dependent AUC for the Joint Model jointFit
#> 
#> Estimated AUC:  0.8374
#> At time: 7
#> Using information up to time: 4 (225 subjects still at risk)
#> Accounting for censoring using model-based weights
#> 
plot(roc, legend = TRUE, optimal_cutoff = "Youden")

# }
```

# JMbayes2 0.5.0

## Major
* `jm()` now allows for zero-correlations constraints in the covariance matrix of the random effects. When the mixed models provided in the `Mixed_objects` argument have been fitted assuming a diagonal matrix for the random effects, this will also be assumed in the joint model (in previous versions, this was ignored). In addition, the new argument `which_independent` can be used to specify which longitudinal outcomes are to be assumed independent.

# JMbayes2 0.4.1

## Major
* Portable implementation of parallel computing.

* function `area()` has gained the argument `time_window` that specifies the window of integrating the linear predictor of the corresponding longitudinal outcome. 


# JMbayes2 0.4.0

## Major
* Function `tvBrier()` has gained the argument `integrated` for calculating the integrated Brier score.

* Function `tvBrier()` has gained the argument `type_weights` and now also allows to correct for censoring in the interval `Tstart` to `Thoriz` using inverse probability of censoring weighting. The default remains model-based weights.

* The new function `tvEPCE()` calculates the time-varying expected predictive cross-entropy.

* This version supports Super Learning for optimizing predictions using cross-validation and a library of joint models. In that regard, the new function `create_folds()` can be used to split a dataset in V-folds of training and test datasets. More information can be found in the corresponding vignette.

## Minor
* Weak informative priors are now used for the fixed-effects of the mixed-effects models.

* Several improvements in various internal functions.

# JMbayes2 0.3.0

## Major
* An issue resulting in wider than expected credible intervals for the fixed-effects coefficients of the longitudinal submodels has been resolved.

## Minor
* Several improvements in various internal functions.

# JMbayes2 0.2.9

## Major
* The default placing of the knots for the B-spline approximation of the log baseline hazard has been changed. This will cause some difference compared to previous versions.

# JMbayes2 0.2.0

## Major
* Dynamic predictions for competing risks data can now be computed. An example is given in the Competing Risks vignette.

* Function `jm()` can now fit joint models with a recurrent event process with or without a terminating event. The model accommodates discontinuous risk intervals, and the time can be defined in terms of the gap or calendar timescale. An example is given in the Recurrent Events vignette.

# JMbayes2 0.1.7

## Major
* Added the function `tvBrier()` for calculating time-varying Brier score for fitted joint models. Currently, only right-censored data are supported.

* Added the functions `calibration_plot()` and `calibration_metrics()` for calculating time-varying calibration plot and calibration metrics for fitted joint models. Currently, only right-censored data are supported.

* Added new section in the vignette for Dynamic Prediction (available on the [website of the package](https://drizopoulos.github.io/JMbayes2/)) to showcase the use of the functions mentioned above. 

## Minor
* Improved the plot method for dynamic predictions.

* Several bug corrections.

# JMbayes2 0.1.6

## Major
* Added a `predict()` method for `jm` objects and a corresponding `plot()` for objects of class `predict_jm` for calculating and displaying predictions from joint models. Currently, only standard survival models are covered. Future versions will include predictions from competing risks and multi-state models.

* Added the functions `tvROC()` and `tvAUC()` for calculating time-varying Receiver Operating Characteristic (ROC) curves and the areas under the ROC curves for fitted joint models. Currently, only right-censored data are supported.

* Added a vignette (available on the [website of the package](https://drizopoulos.github.io/JMbayes2/)) to explain how (dynamic) predictions are calculated in the package. 


# JMbayes2 0.1.5

## Major
* Added two vignettes (available on the [website of the package](https://drizopoulos.github.io/JMbayes2/)) to showcase joint models with competing risks and joint models with non-Gaussian longitudinal outcomes.

* Simplified syntax and additional options for specifying [transformation functions of functional forms](https://drizopoulos.github.io/JMbayes2/articles/Transformation_Functions.html).

* The `slope()` function has gained two new arguments, `eps` and `direction`. This allows calculating the difference of the longitudinal profile over a user-specified interval.


# JMbayes2 0.1.3

## Minor
* Used `parallel::clusterSetRNGStream()` in `jm_fit()` for distributing the seed in the workers.
* Changed the default position of the knots for the B-spline approximation of the log baseline hazard.


# JMbayes2 0.1.2

## Minor
* Changed calls to `floor()` in the C++ code.


# JMbayes2 0.1.0

## General
* First version of the package.


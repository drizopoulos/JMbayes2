# Joint Models for Longitudinal and Time-to-Event Data

Fits multivariate joint models for longitudinal and time-to-event data.

## Usage

``` r
jm(Surv_object, Mixed_objects, time_var, recurrent = FALSE,
  functional_forms = NULL, which_independent = NULL,
  base_hazard = NULL, data_Surv = NULL, id_var = NULL,
  priors = NULL, control = NULL, ...)

value(x, IE_time = NULL)
coefs(x, zero_ind = NULL, IE_time = NULL)
slope(x, eps = 0.001, direction = "both", IE_time = NULL)
velocity(x, eps = 0.001, direction = "both", IE_time = NULL)
acceleration(x, IE_time = NULL)
area(x, time_window = NULL, IE_time = NULL)
Delta(x, time_window = NULL, standardise = TRUE, IE_time = NULL)

vexpit(x)
Dexpit(x)

vexp(x)
Dexp(x)

vabs(x)

vlog(x)
vlog2(x)
vlog10(x)

vsqrt(x)
poly2(x)
poly3(x)
poly4(x)

tv(x, knots = NULL, ord = 2L)
```

## Arguments

- Surv_object:

  an object:\

  - of class 'coxph' fitted by function
    [`coxph()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
    from package **survival**, or

  - of class 'survreg' fitted by function `survreg()` from package
    **survival**.

- Mixed_objects:

  a `list` of objects or a single object. Objects may be:\

  - of class 'lme' fitted by function
    [`lme()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
    from package **nlme**, or

  - of class 'MixMod' fitted by function
    [`mixed_model()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
    from package **GLMMadaptive**.

- time_var:

  a `character` string indicating the time variable in the mixed-effects
  model(s).

- recurrent:

  a `character` string indicating "calendar" or "gap" timescale to fit a
  recurrent event model.

- functional_forms:

  a `list` of formulas. Each formula corresponds to one longitudinal
  outcome and specifies the association structure between that outcome
  and the survival submodel as well as any interaction terms between the
  components of the longitudinal outcome and the survival submodel. See
  **Examples**.

- which_independent:

  a numeric indicator matrix denoting which outcomes are independent. It
  can also be the character string `"all"` in which case all
  longitudinal outcomes are assumed independent. Only relevant in joint
  models with multiple longitudinal outcomes.

- base_hazard:

  a `character` vector indicating the type of hazard function.

- data_Surv:

  the `data.frame` used to fit the Cox/AFT survival submodel.

- id_var:

  a `character` string indicating the id variable in the survival
  submodel.

- priors:

  a named `list` of user-specified prior parameters:

  `mean_betas_HC`

  : the prior mean vector of the normal prior for the regression
    coefficients of the covariates of the longitudinal model(s), which
    were hierarchically centered.

  `Tau_betas_HC`

  : the prior precision matrix of the normal prior for the regression
    coefficients of the longitudinal model(s), which were hierarchically
    centered.

  `mean_betas_nHC`

  : a `list` of the prior mean vector(s) of the normal prior(s) for the
    regression coefficients of the covariates of the longitudinal
    model(s), which were not hierarchically centered.

  `Tau_betas_nHC`

  : a `list` of the prior precision matrix(ces) of the normal prior(s)
    for the regression coefficients of the longitudinal model(s), which
    were not Hierarchically Centered.

  `mean_bs_gammas`

  : the prior mean vector of the normal prior for the B-splines
    coefficients used to approximate the baseline hazard.

  `Tau_bs_gammas`

  : the prior precision matrix of the normal prior for the B-splines
    coefficients used to approximate the baseline hazard.

  `A_tau_bs_gammas`

  : the prior shape parameter of the gamma prior for the precision
    parameter of the penalty term for the B-splines coefficients for the
    baseline hazard.

  `B_tau_bs_gammas`

  : the prior rate parameter of the gamma prior for the precision
    parameter of the penalty term for the B-splines coefficients for the
    baseline hazard.

  `rank_Tau_bs_gammas`

  : the prior rank parameter for the precision matrix of the normal
    prior for the B-splines coefficients used to approximate the
    baseline hazard.

  `mean_gammas`

  : the prior mean vector of the normal prior for the regression
    coefficients of baseline covariates.

  `Tau_gammas`

  : the prior precision matrix of the normal prior for the regression
    coefficients of baseline covariates.

  `penalty_gammas`

  : a character string with value 'none', 'ridge', or 'horseshoe'
    indicating whether the coefficients of the baseline covariates
    included in the survival submodel should not be shrunk, shrank using
    ridge prior, or shrank using horseshoe prior, respectively.

  `A_lambda_gammas`

  : the prior shape parameter of the gamma prior for the precision
    parameter of the local penalty term for the baseline regression
    coefficients. Only relevant when `penalty_gammas = 'ridge'` or when
    `penalty_gammas = 'horseshoe'`.

  `B_lambda_gammas`

  : the prior rate parameter of the gamma prior for the precision
    parameter of the local penalty term for the baseline regression
    coefficients. Only relevant when `penalty_gammas = 'ridge'` or when
    `penalty_gammas = 'horseshoe'`.

  `A_tau_gammas`

  : the prior shape parameter of the gamma prior for the precision
    parameter of the global penalty term for the baseline regression
    coefficients. Only relevant when `penalty_gammas = 'ridge'` or when
    `penalty_gammas = 'horseshoe'`.

  `B_tau_gammas`

  : the prior rate parameter of the gamma prior for the precision
    parameter of the global penalty term for the baseline regression
    coefficients. Only relevant when `penalty_gammas = 'ridge'` or when
    `penalty_gammas = 'horseshoe'`.

  `A_nu_gammas`

  : the prior shape parameter of the gamma prior for the variance
    hyperparameter for the precision parameter of the local penalty term
    for the baseline regression coefficients. Only relevant when
    `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `B_nu_gammas`

  : the prior rate parameter of the gamma prior for the variance
    hyperparameter for the precision parameter of the local penalty term
    for the baseline regression coefficients. Only relevant when
    `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `A_xi_gammas`

  : the prior shape parameter of the gamma prior for the variance
    hyperparameter for the precision parameter of the global penalty
    term for the baseline regression coefficients. Only relevant when
    `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `B_xi_gammas`

  : the prior rate parameter of the gamma prior for the variance
    hyperparameter for the precision parameter of the global penalty
    term for the baseline regression coefficients. Only relevant when
    `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `mean_alphas`

  : the prior mean vector of the normal prior for the association
    parameter(s).

  `Tau_alphas`

  : the prior mean vector of the normal prior for the association
    parameter(s).

  `penalty_alphas`

  : a character string with value 'none', 'ridge', 'horseshoe'
    indicating whether the coefficients association parameters should
    not be shrunk, shrank using ridge prior, or shrank using horseshoe
    prior, respectively.

  `A_lambda_alphas`

  : the prior shape parameter of the gamma prior for the precision
    parameter of the local penalty term for the association parameters.
    Only relevant when `penalty_gammas = 'ridge'` or when
    `penalty_gammas = 'horseshoe'`.

  `B_lambda_alphas`

  : the prior rate parameter of the gamma prior for the precision
    parameter of the local penalty term for the association parameters.
    Only relevant when `penalty_gammas = 'ridge'` or when
    `penalty_gammas = 'horseshoe'`.

  `A_tau_alphas`

  : the prior shape parameter of the gamma prior for the precision
    parameter of the global penalty term for the association parameters.
    Only relevant when `penalty_gammas = 'ridge'` or when
    `penalty_gammas = 'horseshoe'`.

  `B_tau_alphas`

  : the prior rate parameter of the gamma prior for the precision
    parameter of the global penalty term for the association parameters.
    Only relevant when `penalty_gammas = 'ridge'` or
    `penalty_gammas = 'horseshoe'`.

  `A_nu_alphas`

  : the prior shape parameter of the gamma prior for the variance
    hyperparameter for the precision parameter of the local penalty term
    for the association parameters. Only relevant when
    `penalty_gammas = 'ridge'`, or `penalty_gammas = 'horseshoe'`.

  `B_nu_alphas`

  : the prior rate parameter of the gamma prior for the variance
    hyperparameter for the precision parameter of the local penalty term
    for the association parameters. Only relevant when
    `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `A_xi_alphas`

  : the prior shape parameter of the gamma prior for the variance
    hyperparameter for the precision parameter of the global penalty
    term for the association parameters. Only relevant when
    `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `B_xi_alphas`

  : the prior rate parameter of the gamma prior for the variance
    hyperparameter for the precision parameter of the global penalty
    term for the association parameters. Only relevant when
    `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `gamma_prior_D_sds`

  : logical; if `TRUE`, a gamma prior will be used for the standard
    deviations of the D matrix (variance-covariance matrix of the random
    effects). Defaults to `TRUE`

  `D_sds_df`

  : the prior degrees of freedom parameter for the half-t prior for the
    standard deviations of the D matrix (variance-covariance matrix of
    the random effects).

  `D_sds_sigma`

  : the prior sigma parameter vector for the half-t prior for the
    standard deviations of the D matrix (variance-covariance matrix of
    the random effects).

  `D_sds_shape`

  : the prior shape parameter for the gamma prior for the standard
    deviations of the D matrix (variance-covariance matrix of the random
    effects).

  `D_sds_mean`

  : the prior mean parameter vector for the gamma prior for the standard
    deviations of the D matrix (variance-covariance matrix of the random
    effects).

  `D_L_etaLKJ`

  : the prior eta parameter for the LKJ prior for the correlation matrix
    of the random effects.

  `sigmas_df`

  : the prior degrees of freedom parameter for the half-t prior for the
    error term(s).

  `sigmas_sigma`

  : the prior sigma parameter for the half-t prior for the error
    term(s).

- control:

  a list of control values with components:

  `GK_k`

  : the number of quadrature points for the Gauss Kronrod rule; options
    15 and 7.

  `n_chains`

  : an integer specifying the number of chains for the MCMC. Defaults to
    3.

  `n_burnin`

  : an integer specifying the number of burn-in iterations. Defaults to
    500.

  `n_iter`

  : an integer specifying the number of total iterations per chain.
    Defaults to 3500.

  `n_thin`

  : an integer specifying the thinning of the chains. Defaults to 1.

  `seed`

  : the seed used in the sampling procedures. Defaults to 123.

  `MALA`

  : `logical`; if `TRUE`, the MALA algorithm is used when updating the
    elements of the Cholesky factor of the D matrix. Defaults to
    `FALSE`.

  `save_random_effects`

  : `logical`; if `TRUE`, the full MCMC results of the random effects
    will be saved and returned with the `jm` object. Defaults to
    `FALSE`.

  `save_logLik_contributions`

  : `logical`; if `TRUE`, the log-likelihood contributions are saved in
    the `mcmc` component of the `jm` object. Defaults to `FALSE`

  `cores`

  : an integer specifying the number of cores to use for running the
    chains in parallel; no point of setting this greater than
    `n_chains`.

  `parallel`

  : a character string indicating how the parallel sampling of the
    chains will be performed. Options are `"snow"` (default) and
    `"multicore"`.

  `basis`

  : character string with possible values `"bs"` (default) or `"ns"`.
    When `"bs"` a B-spline basis is used to approximate the log baseline
    hazard function with degree of the spline specified by the
    `Bsplines_degree`. When `"ns"` a natrual cubic spline basis is used;
    in this case the value of the `Bsplines_degree` control argument is
    ignored.

  `Bsplines_degree`

  : the degree of the splines in each basis; default is quadratic
    splines.

  `base_hazard_segments`

  : the number of segments to split the follow-up period for the spline
    approximation of the log baseline hazard function. Defaults to 10.

  `timescale_base_hazard`

  : character string with possible values `"identity"` (default) or
    `"log"`. When `"identity"` the spline basis is specified for the
    time variable in its orginal scale. When `"log"` the spline basis is
    specified for the logarithm of the time variable.

  `diff`

  : the order of the difference used in the penalty matrix for the
    coefficients of the splines used to approximate the log baseline
    hazard function. Defaults to 2.

  `knots`

  : a numeric vector with the position of the knots for the spline
    approximation of the log baseline hazard function. The default is
    equally-spaced knots starting from `sqrt(.Machine$double.eps)` until
    the maximum follow-up time.

- x:

  a numeric input variable.

- knots:

  a numeric vector of knots.

- ord:

  an integer denoting the order of the spline.

- zero_ind:

  a list with integer vectors indicating which coefficients are set to
  zero in the calculation of the value term. This can be used to include
  for example only the random intercept; default is `NULL`.

- eps:

  numeric scalar denoting the step-size for the finite difference
  approximation.

- direction:

  character string for the direction of the numerical derivative,
  options are `"both"`, and `"backward"`.

- time_window:

  numeric scalar controlling the time window used by the `Delta()` and
  `area()` functional forms. For `area()`, `time_window` specifies the
  lower limit of the interval over which the integral is evaluated. For
  `Delta()`, `time_window` specifies the length of the time interval
  (i.e., the contrast between t and t - `time_window`) over which the
  finite difference is computed; when set to `NULL` (the default), the
  contrast is taken between the current time and time 0.

- standardise:

  `logical`; controls whether the `Delta()` functional form returns a
  rate or a raw contrast. If `TRUE`, the difference between the values
  at the two time points is divided by the time distance between them,
  yielding a change per unit time. If `FALSE`, `Delta()` returns the raw
  difference over the specified time window. Defaults to `TRUE`.

- IE_time:

  a `character` string specifying the name of the intermediate event
  time variable in the `data.frame` used to fit the Cox/AFT survival
  submodel. For groups/subjects who did not experience the intermediate
  event, the time should be set to `Inf`. The same `IE_time` variable
  should be used when specifying multiple functional forms for the same
  longitudinal outcome.

- ...:

  arguments passed to `control`.

## Details

The mathematical details regarding the definition of the multivariate
joint model, and the capabilities of the package can be found in the
vignette in the doc directory.

Notes:

- The ordering of the subjects in the datasets used to fit the mixed and
  Cox regression models needs to be the same.

- The units of the time variables in the mixed and Cox models need to be
  the same.

## Value

A list of class `jm` with components:

- mcmc:

  a `list` of the MCMC samples for each parameter.

- acc_rates:

  a `list` of the acceptance rates for each parameter.

- logLik:

  a `matrix` of dimensions \[`((n_iter - n_burnin)/n_thin)*n_thin`,
  number of individuals\], with element \[i, j\] being the conditional
  log-Likelihood value of the \\i^{th}\\ iteration for the \\j^{th}\\
  individual.

- mlogLik:

  a `matrix` of dimensions \[`((n_iter - n_burnin)/n_thin)*n_thin`,
  number of individuals\], with element \[i, j\] being the marginal
  log-Likelihood value of the \\i^{th}\\ iteration for the \\j^{th}\\
  individual.

- running_time:

  an object of class `proc_time` with the time used to run `jm`.

- statistics:

  a `list` with posterior estimates of the parameters (means, medians,
  standard deviations, standard errors, effective sample sizes, tail
  probabilities, upper and lower bounds of credible intervals, etc.).

- fit_stats:

  a `list` of lists with fit statistics (DIC, pD, LPML, CPO, WAIC) for
  both conditional and marginal formulations.

- model_data:

  a `list` of data used to fit the model.

- model_info:

  a `list` of components of the fit useful to other functions.

- initial_values:

  a `list` with the initial values of the parameters.

- control:

  a copy of the `control` values used to fit the model.

- priors:

  a copy of the `priors` used to fit the model.

- call:

  the matched call.

## Author

Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>

## See also

[`methods.jm`](https://drizopoulos.github.io/JMbayes2/reference/methods.md),
[`coda_methods.jm`](https://drizopoulos.github.io/JMbayes2/reference/coda_methods.md)

## Examples

``` r
# \donttest{
################################################################################

##############################################
# Univariate joint model for serum bilirubin #
# 1 continuous outcome                       #
##############################################

# [1] Fit the mixed model using lme().
fm1 <- lme(fixed = log(serBilir) ~ year * sex + I(year^2) +
           age + prothrombin, random =  ~ year | id, data = pbc2)

# [2] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [3] The basic joint model is fitted using a call to jm() i.e.,
joint_model_fit_1 <- jm(fCox1, fm1, time_var = "year",
        n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)
summary(joint_model_fit_1)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = fCox1, Mixed_objects = fm1, time_var = "year", 
#>     n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)
#> 
#> Data Descriptives:
#> Number of groups: 312        Number of events: 140 (44.9%)
#> Number of observations:
#>   log(serBilir): 1945
#> 
#>                  DIC     WAIC      LPML
#> marginal    4205.481 5053.292 -3002.078
#> conditional 3341.875 3171.322 -1806.821
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.9782 (Intr)
#> year   0.1769 0.3433
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%      P
#> drugD-penicil        -0.0250 0.2362 -0.4807 0.4442 0.9112
#> age                   0.0638 0.0089  0.0461 0.0814 0.0000
#> value(log(serBilir))  1.4240 0.1058  1.2261 1.6298 0.0000
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P
#> (Intercept)     0.2542 0.3571 -0.4535  0.9539 0.4754
#> year            0.2265 0.0371  0.1546  0.2999 0.0000
#> sexfemale      -0.2448 0.1805 -0.5969  0.1096 0.1718
#> I(year^2)       0.0026 0.0011  0.0003  0.0047 0.0262
#> age            -0.0018 0.0054 -0.0123  0.0088 0.7340
#> prothrombin     0.0527 0.0085  0.0364  0.0690 0.0000
#> year:sexfemale -0.0866 0.0385 -0.1628 -0.0131 0.0240
#> sigma           0.3456 0.0068  0.3326  0.3594 0.0000
#> 
#> MCMC summary:
#> chains: 1 
#> iterations per chain: 11000 
#> burn-in per chain: 1000 
#> thinning: 1 
#> time: 23 sec
traceplot(joint_model_fit_1)



























################################################################################

##########################################################################
# Multivariate joint model for serum bilirubin, hepatomegaly and ascites #
# 1 continuous outcome, 2 categorical outcomes                           #
##########################################################################

# [1] Fit the mixed-effects models using lme() for continuous
# outcomes and mixed_model() for categorical outcomes.
fm1 <- lme(fixed = log(serBilir) ~ year * sex,
           random = ~ year | id, data = pbc2)

fm2 <- mixed_model(hepatomegaly ~ sex + age + year, data = pbc2,
                   random = ~ year | id, family = binomial())

fm3 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ year | id, family = binomial())

# [2] Save all the fitted mixed-effects models in a list.
Mixed <- list(fm1, fm2, fm3)

# [3] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [4] The joint model is fitted using a call to jm() i.e.,
joint_model_fit_2 <- jm(fCox1, Mixed, time_var = "year",
      n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)
summary(joint_model_fit_2)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = fCox1, Mixed_objects = Mixed, time_var = "year", 
#>     n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)
#> 
#> Data Descriptives:
#> Number of groups: 312        Number of events: 140 (44.9%)
#> Number of observations:
#>   log(serBilir): 1945
#>   hepatomegaly: 1884
#>   ascites: 1885
#> 
#>                  DIC     WAIC      LPML
#> marginal    6637.755 6900.663 -3849.614
#> conditional 9091.606 8804.787 -4827.846
#> 
#> Random-effects covariance matrix:
#>                                                   
#>        StdDev   Corr                              
#> (Intr) 0.9907 (Intr)   year (Intr)    year  (Intr)
#> year   0.1780 0.4053                              
#> (Intr) 3.3854 0.5323 0.3526                       
#> year   0.5905 0.0457 0.3445 -0.3450               
#> (Intr) 2.7652 0.6224 0.5219 0.5671  -0.0528       
#> year   0.4153 0.4182 0.6459 0.3522  0.3312  0.0235
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%      P
#> drugD-penicil        -0.2043 0.2763 -0.7681 0.3237 0.4550
#> age                   0.0334 0.0148  0.0022 0.0601 0.0364
#> value(log(serBilir))  0.6809 0.2248  0.2445 1.1019 0.0000
#> value(hepatomegaly)  -0.0754 0.0869 -0.2535 0.0954 0.3820
#> value(ascites)        0.6089 0.1943  0.2163 0.9939 0.0002
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P
#> (Intercept)     0.6614 0.1402  0.3879  0.9396 0.0000
#> year            0.2493 0.0288  0.1940  0.3071 0.0000
#> sexfemale      -0.1935 0.1428 -0.4746  0.0836 0.1702
#> year:sexfemale -0.0645 0.0297 -0.1246 -0.0078 0.0274
#> sigma           0.3478 0.0067  0.3351  0.3615 0.0000
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>                Mean  StDev    2.5%  97.5%      P
#> (Intercept)  0.0919 1.0555 -1.9639 2.1926 0.9378
#> sexfemale   -0.7664 0.5335 -1.8295 0.2783 0.1424
#> age          0.0147 0.0170 -0.0185 0.0483 0.3834
#> year         0.2537 0.0757  0.1147 0.4182 0.0000
#> 
#> Longitudinal outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev     2.5%   97.5% P
#> (Intercept) -8.8247 0.9669 -10.8540 -7.0448 0
#> year         0.5704 0.0587   0.4546  0.6856 0
#> age          0.0796 0.0162   0.0502  0.1124 0
#> 
#> MCMC summary:
#> chains: 1 
#> iterations per chain: 11000 
#> burn-in per chain: 1000 
#> thinning: 1 
#> time: 1 min
traceplot(joint_model_fit_2)



















































################################################################################

######################
# Slope & Area Terms #
######################

# We extend model 'joint_model_fit_2' by including the value and slope term for
# bilirubin, the area term for hepatomegaly (in the log-odds scale), and the
# value and area term for spiders (in the log-odds scale).
# To include these terms into the model, we specify the 'functional_forms'
# argument. This should be a list of right side formulas. Each component of the
# list should have as name the name of the corresponding outcome variable. In
# the right side formula we specify the functional form of the association using
# functions 'value()', 'slope()' and 'area()'.
# Notes: (1) For terms not specified in the 'functional_forms' list, the default
# value functional form is used.

# [1] Fit the mixed-effects models using lme() for continuous outcomes
# and mixed_model() for categorical outcomes.
fm1 <- lme(fixed = log(serBilir) ~ year * sex, random = ~ year | id, data = pbc2)

fm2 <- mixed_model(hepatomegaly ~ sex + age + year, data = pbc2,
                   random = ~ year | id, family = binomial())

fm3 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ year | id, family = binomial())

# [2] Save all the fitted mixed-effects models in a list.
Mixed <- list(fm1, fm2, fm3)

# [3] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [4] Specify the list of formulas to be passed to the functional_forms argument
# of jm().
fForms <- list("log(serBilir)" = ~ value(log(serBilir)) + slope(log(serBilir)),
               "hepatomegaly" = ~ area(hepatomegaly),
               "ascites" = ~ value(ascites) + area(ascites))

# [5] The joint model is fitted using a call to jm() and passing the list
# to the functional_forms argument.
joint_model_fit_2 <- jm(fCox1, Mixed, time_var = "year",
                        functional_forms = fForms, n_chains = 1L,
                        n_iter = 11000L, n_burnin = 1000L)
summary(joint_model_fit_2)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = fCox1, Mixed_objects = Mixed, time_var = "year", 
#>     functional_forms = fForms, n_chains = 1L, n_iter = 11000L, 
#>     n_burnin = 1000L)
#> 
#> Data Descriptives:
#> Number of groups: 312        Number of events: 140 (44.9%)
#> Number of observations:
#>   log(serBilir): 1945
#>   hepatomegaly: 1884
#>   ascites: 1885
#> 
#>                  DIC     WAIC     LPML
#> marginal    6649.297 6979.906 -3987.59
#> conditional 9027.793 8783.470 -4853.26
#> 
#> Random-effects covariance matrix:
#>                                                   
#>        StdDev   Corr                              
#> (Intr) 0.9838 (Intr)   year (Intr)    year  (Intr)
#> year   0.1800 0.4185                              
#> (Intr) 3.3705 0.5314 0.3362                       
#> year   0.6092 0.0524 0.3824 -0.3449               
#> (Intr) 2.3196 0.6533 0.5711 0.5921  -0.0850       
#> year   0.4023 0.4636 0.6857 0.4244  0.3691  0.3019
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%   97.5%      P
#> drugD-penicil        -0.1273 0.2727 -0.6748  0.3839 0.6372
#> age                   0.0476 0.0133  0.0210  0.0730 0.0008
#> value(log(serBilir))  0.8984 0.1934  0.5395  1.2775 0.0000
#> slope(log(serBilir))  3.9438 1.3701  1.1925  6.8787 0.0018
#> area(hepatomegaly)    0.1276 0.0749 -0.0241  0.2838 0.0868
#> value(ascites)       -0.6944 0.1906 -1.0284 -0.2922 0.0000
#> area(ascites)         1.0414 0.2576  0.5694  1.5365 0.0000
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P
#> (Intercept)     0.7087 0.1439  0.4292  0.9893 0.0000
#> year            0.2599 0.0291  0.2052  0.3190 0.0000
#> sexfemale      -0.2472 0.1480 -0.5448  0.0411 0.0904
#> year:sexfemale -0.0716 0.0297 -0.1310 -0.0144 0.0128
#> sigma           0.3483 0.0066  0.3357  0.3615 0.0000
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>                Mean  StDev    2.5%  97.5%      P
#> (Intercept)  0.2779 1.0515 -1.7473 2.3695 0.8068
#> sexfemale   -0.9298 0.5511 -2.0285 0.1564 0.0920
#> age          0.0138 0.0167 -0.0196 0.0459 0.4034
#> year         0.2688 0.0705  0.1388 0.4133 0.0002
#> 
#> Longitudinal outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev     2.5%   97.5% P
#> (Intercept) -8.2935 0.9795 -10.2808 -6.4226 0
#> year         0.4666 0.0721   0.3209  0.6113 0
#> age          0.0762 0.0164   0.0447  0.1090 0
#> 
#> MCMC summary:
#> chains: 1 
#> iterations per chain: 11000 
#> burn-in per chain: 1000 
#> thinning: 1 
#> time: 1.1 min

# }
```

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

  an object:  

  - of class 'coxph' fitted by function
    [`coxph()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
    from package **survival**, or

  - of class 'survreg' fitted by function `survreg()` from package
    **survival**.

- Mixed_objects:

  a `list` of objects or a single object. Objects may be:  

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

  :   the prior mean vector of the normal prior for the regression
      coefficients of the covariates of the longitudinal model(s), which
      were hierarchically centered.

  `Tau_betas_HC`

  :   the prior precision matrix of the normal prior for the regression
      coefficients of the longitudinal model(s), which were
      hierarchically centered.

  `mean_betas_nHC`

  :   a `list` of the prior mean vector(s) of the normal prior(s) for
      the regression coefficients of the covariates of the longitudinal
      model(s), which were not hierarchically centered.

  `Tau_betas_nHC`

  :   a `list` of the prior precision matrix(ces) of the normal prior(s)
      for the regression coefficients of the longitudinal model(s),
      which were not Hierarchically Centered.

  `mean_bs_gammas`

  :   the prior mean vector of the normal prior for the B-splines
      coefficients used to approximate the baseline hazard.

  `Tau_bs_gammas`

  :   the prior precision matrix of the normal prior for the B-splines
      coefficients used to approximate the baseline hazard.

  `A_tau_bs_gammas`

  :   the prior shape parameter of the gamma prior for the precision
      parameter of the penalty term for the B-splines coefficients for
      the baseline hazard.

  `B_tau_bs_gammas`

  :   the prior rate parameter of the gamma prior for the precision
      parameter of the penalty term for the B-splines coefficients for
      the baseline hazard.

  `rank_Tau_bs_gammas`

  :   the prior rank parameter for the precision matrix of the normal
      prior for the B-splines coefficients used to approximate the
      baseline hazard.

  `mean_gammas`

  :   the prior mean vector of the normal prior for the regression
      coefficients of baseline covariates.

  `Tau_gammas`

  :   the prior precision matrix of the normal prior for the regression
      coefficients of baseline covariates.

  `penalty_gammas`

  :   a character string with value 'none', 'ridge', or 'horseshoe'
      indicating whether the coefficients of the baseline covariates
      included in the survival submodel should not be shrunk, shrank
      using ridge prior, or shrank using horseshoe prior, respectively.

  `A_lambda_gammas`

  :   the prior shape parameter of the gamma prior for the precision
      parameter of the local penalty term for the baseline regression
      coefficients. Only relevant when `penalty_gammas = 'ridge'` or
      when `penalty_gammas = 'horseshoe'`.

  `B_lambda_gammas`

  :   the prior rate parameter of the gamma prior for the precision
      parameter of the local penalty term for the baseline regression
      coefficients. Only relevant when `penalty_gammas = 'ridge'` or
      when `penalty_gammas = 'horseshoe'`.

  `A_tau_gammas`

  :   the prior shape parameter of the gamma prior for the precision
      parameter of the global penalty term for the baseline regression
      coefficients. Only relevant when `penalty_gammas = 'ridge'` or
      when `penalty_gammas = 'horseshoe'`.

  `B_tau_gammas`

  :   the prior rate parameter of the gamma prior for the precision
      parameter of the global penalty term for the baseline regression
      coefficients. Only relevant when `penalty_gammas = 'ridge'` or
      when `penalty_gammas = 'horseshoe'`.

  `A_nu_gammas`

  :   the prior shape parameter of the gamma prior for the variance
      hyperparameter for the precision parameter of the local penalty
      term for the baseline regression coefficients. Only relevant when
      `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `B_nu_gammas`

  :   the prior rate parameter of the gamma prior for the variance
      hyperparameter for the precision parameter of the local penalty
      term for the baseline regression coefficients. Only relevant when
      `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `A_xi_gammas`

  :   the prior shape parameter of the gamma prior for the variance
      hyperparameter for the precision parameter of the global penalty
      term for the baseline regression coefficients. Only relevant when
      `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `B_xi_gammas`

  :   the prior rate parameter of the gamma prior for the variance
      hyperparameter for the precision parameter of the global penalty
      term for the baseline regression coefficients. Only relevant when
      `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `mean_alphas`

  :   the prior mean vector of the normal prior for the association
      parameter(s).

  `Tau_alphas`

  :   the prior mean vector of the normal prior for the association
      parameter(s).

  `penalty_alphas`

  :   a character string with value 'none', 'ridge', 'horseshoe'
      indicating whether the coefficients association parameters should
      not be shrunk, shrank using ridge prior, or shrank using horseshoe
      prior, respectively.

  `A_lambda_alphas`

  :   the prior shape parameter of the gamma prior for the precision
      parameter of the local penalty term for the association
      parameters. Only relevant when `penalty_gammas = 'ridge'` or when
      `penalty_gammas = 'horseshoe'`.

  `B_lambda_alphas`

  :   the prior rate parameter of the gamma prior for the precision
      parameter of the local penalty term for the association
      parameters. Only relevant when `penalty_gammas = 'ridge'` or when
      `penalty_gammas = 'horseshoe'`.

  `A_tau_alphas`

  :   the prior shape parameter of the gamma prior for the precision
      parameter of the global penalty term for the association
      parameters. Only relevant when `penalty_gammas = 'ridge'` or when
      `penalty_gammas = 'horseshoe'`.

  `B_tau_alphas`

  :   the prior rate parameter of the gamma prior for the precision
      parameter of the global penalty term for the association
      parameters. Only relevant when `penalty_gammas = 'ridge'` or
      `penalty_gammas = 'horseshoe'`.

  `A_nu_alphas`

  :   the prior shape parameter of the gamma prior for the variance
      hyperparameter for the precision parameter of the local penalty
      term for the association parameters. Only relevant when
      `penalty_gammas = 'ridge'`, or `penalty_gammas = 'horseshoe'`.

  `B_nu_alphas`

  :   the prior rate parameter of the gamma prior for the variance
      hyperparameter for the precision parameter of the local penalty
      term for the association parameters. Only relevant when
      `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `A_xi_alphas`

  :   the prior shape parameter of the gamma prior for the variance
      hyperparameter for the precision parameter of the global penalty
      term for the association parameters. Only relevant when
      `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `B_xi_alphas`

  :   the prior rate parameter of the gamma prior for the variance
      hyperparameter for the precision parameter of the global penalty
      term for the association parameters. Only relevant when
      `penalty_gammas = 'ridge'` or when `penalty_gammas = 'horseshoe'`.

  `gamma_prior_D_sds`

  :   logical; if `TRUE`, a gamma prior will be used for the standard
      deviations of the D matrix (variance-covariance matrix of the
      random effects). Defaults to `TRUE`

  `D_sds_df`

  :   the prior degrees of freedom parameter for the half-t prior for
      the standard deviations of the D matrix (variance-covariance
      matrix of the random effects).

  `D_sds_sigma`

  :   the prior sigma parameter vector for the half-t prior for the
      standard deviations of the D matrix (variance-covariance matrix of
      the random effects).

  `D_sds_shape`

  :   the prior shape parameter for the gamma prior for the standard
      deviations of the D matrix (variance-covariance matrix of the
      random effects).

  `D_sds_mean`

  :   the prior mean parameter vector for the gamma prior for the
      standard deviations of the D matrix (variance-covariance matrix of
      the random effects).

  `D_L_etaLKJ`

  :   the prior eta parameter for the LKJ prior for the correlation
      matrix of the random effects.

  `sigmas_df`

  :   the prior degrees of freedom parameter for the half-t prior for
      the error term(s).

  `sigmas_sigma`

  :   the prior sigma parameter for the half-t prior for the error
      term(s).

- control:

  a list of control values with components:

  `GK_k`

  :   the number of quadrature points for the Gauss Kronrod rule;
      options 15 and 7.

  `n_chains`

  :   an integer specifying the number of chains for the MCMC. Defaults
      to 3.

  `n_burnin`

  :   an integer specifying the number of burn-in iterations. Defaults
      to 500.

  `n_iter`

  :   an integer specifying the number of total iterations per chain.
      Defaults to 3500.

  `n_thin`

  :   an integer specifying the thinning of the chains. Defaults to 1.

  `seed`

  :   the seed used in the sampling procedures. Defaults to 123.

  `MALA`

  :   `logical`; if `TRUE`, the MALA algorithm is used when updating the
      elements of the Cholesky factor of the D matrix. Defaults to
      `FALSE`.

  `save_random_effects`

  :   `logical`; if `TRUE`, the full MCMC results of the random effects
      will be saved and returned with the `jm` object. Defaults to
      `FALSE`.

  `save_logLik_contributions`

  :   `logical`; if `TRUE`, the log-likelihood contributions are saved
      in the `mcmc` component of the `jm` object. Defaults to `FALSE`

  `cores`

  :   an integer specifying the number of cores to use for running the
      chains in parallel; no point of setting this greater than
      `n_chains`.

  `parallel`

  :   a character string indicating how the parallel sampling of the
      chains will be performed. Options are `"snow"` (default) and
      `"multicore"`.

  `basis`

  :   character string with possible values `"bs"` (default) or `"ns"`.
      When `"bs"` a B-spline basis is used to approximate the log
      baseline hazard function with degree of the spline specified by
      the `Bsplines_degree`. When `"ns"` a natrual cubic spline basis is
      used; in this case the value of the `Bsplines_degree` control
      argument is ignored.

  `Bsplines_degree`

  :   the degree of the splines in each basis; default is quadratic
      splines.

  `base_hazard_segments`

  :   the number of segments to split the follow-up period for the
      spline approximation of the log baseline hazard function. Defaults
      to 10.

  `timescale_base_hazard`

  :   character string with possible values `"identity"` (default) or
      `"log"`. When `"identity"` the spline basis is specified for the
      time variable in its orginal scale. When `"log"` the spline basis
      is specified for the logarithm of the time variable.

  `diff`

  :   the order of the difference used in the penalty matrix for the
      coefficients of the splines used to approximate the log baseline
      hazard function. Defaults to 2.

  `knots`

  :   a numeric vector with the position of the knots for the spline
      approximation of the log baseline hazard function. The default is
      equally-spaced knots starting from `sqrt(.Machine$double.eps)`
      until the maximum follow-up time.

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
#> marginal    4204.643 5075.364 -2931.823
#> conditional 3334.409 3165.829 -1817.475
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.9752 (Intr)
#> year   0.1772 0.3429
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%     P
#> drugD-penicil        -0.0294 0.2359 -0.4760 0.4417 0.913
#> age                   0.0639 0.0092  0.0465 0.0824 0.000
#> value(log(serBilir))  1.4254 0.0976  1.2296 1.6264 0.000
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P
#> (Intercept)     0.2444 0.3605 -0.4688  0.9355 0.4880
#> year            0.2281 0.0370  0.1568  0.3005 0.0000
#> sexfemale      -0.2421 0.1816 -0.6013  0.1134 0.1850
#> I(year^2)       0.0026 0.0010  0.0007  0.0045 0.0066
#> age            -0.0017 0.0054 -0.0122  0.0092 0.7452
#> prothrombin     0.0529 0.0085  0.0364  0.0695 0.0000
#> year:sexfemale -0.0881 0.0385 -0.1644 -0.0137 0.0206
#> sigma           0.3452 0.0068  0.3322  0.3594 0.0000
#> 
#> MCMC summary:
#> chains: 1 
#> iterations per chain: 11000 
#> burn-in per chain: 1000 
#> thinning: 1 
#> time: 29 sec
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
#> marginal    6632.739 6870.898 -3753.886
#> conditional 9110.856 8863.104 -4876.052
#> 
#> Random-effects covariance matrix:
#>                                                    
#>        StdDev   Corr                               
#> (Intr) 0.9893 (Intr)   year (Intr)    year  (Intr) 
#> year   0.1733 0.3907                               
#> (Intr) 3.2864 0.5330 0.3398                        
#> year   0.5734 0.0382 0.3498 -0.3487                
#> (Intr) 2.9721 0.6096 0.4818 0.5282  -0.0083        
#> year   0.4331 0.3691 0.5934 0.3360  0.2505  -0.1014
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%      P
#> drugD-penicil        -0.1835 0.2681 -0.6967 0.3418 0.4940
#> age                   0.0345 0.0138  0.0055 0.0602 0.0228
#> value(log(serBilir))  0.7152 0.2214  0.2396 1.1226 0.0022
#> value(hepatomegaly)  -0.0422 0.0895 -0.2327 0.1400 0.6166
#> value(ascites)        0.5345 0.1802  0.1911 0.9177 0.0000
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P
#> (Intercept)     0.6763 0.1390  0.4064  0.9513 0.0000
#> year            0.2449 0.0274  0.1923  0.2990 0.0000
#> sexfemale      -0.2085 0.1435 -0.4946  0.0675 0.1442
#> year:sexfemale -0.0641 0.0281 -0.1203 -0.0106 0.0182
#> sigma           0.3482 0.0067  0.3353  0.3618 0.0000
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>                Mean  StDev    2.5%  97.5%      P
#> (Intercept)  0.1123 1.0239 -1.8832 2.1349 0.9162
#> sexfemale   -0.8010 0.5274 -1.8554 0.1974 0.1188
#> age          0.0146 0.0164 -0.0176 0.0470 0.3742
#> year         0.2485 0.0694  0.1149 0.3874 0.0000
#> 
#> Longitudinal outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev     2.5%   97.5% P
#> (Intercept) -9.1784 1.1370 -11.4720 -7.1558 0
#> year         0.5916 0.0751   0.4611  0.7529 0
#> age          0.0829 0.0172   0.0507  0.1169 0
#> 
#> MCMC summary:
#> chains: 1 
#> iterations per chain: 11000 
#> burn-in per chain: 1000 
#> thinning: 1 
#> time: 1.1 min
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
#>                  DIC     WAIC      LPML
#> marginal    6687.050 7034.445 -4283.997
#> conditional 9002.434 8739.099 -4827.833
#> 
#> Random-effects covariance matrix:
#>                                                   
#>        StdDev   Corr                              
#> (Intr) 0.9891 (Intr)   year (Intr)    year  (Intr)
#> year   0.1798 0.4155                              
#> (Intr) 3.3535 0.5361 0.3412                       
#> year   0.5831 0.0243 0.3575 -0.3801               
#> (Intr) 2.2908 0.6482 0.5783 0.6036  -0.1252       
#> year   0.4029 0.4939 0.6868 0.4348  0.3518  0.3231
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%   97.5%      P
#> drugD-penicil        -0.1218 0.2734 -0.6696  0.3927 0.6718
#> age                   0.0472 0.0132  0.0188  0.0708 0.0006
#> value(log(serBilir))  0.9030 0.1869  0.5548  1.2815 0.0000
#> slope(log(serBilir))  3.6008 1.2959  1.0806  6.2951 0.0034
#> area(hepatomegaly)    0.1217 0.0867 -0.0535  0.3053 0.1558
#> value(ascites)       -0.6116 0.2099 -0.9474 -0.1067 0.0034
#> area(ascites)         0.9588 0.2968  0.1782  1.4429 0.0028
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P
#> (Intercept)     0.7047 0.1433  0.4263  0.9870 0.0000
#> year            0.2590 0.0294  0.2023  0.3175 0.0000
#> sexfemale      -0.2421 0.1481 -0.5327  0.0457 0.0998
#> year:sexfemale -0.0715 0.0302 -0.1323 -0.0136 0.0152
#> sigma           0.3481 0.0067  0.3355  0.3617 0.0000
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>                Mean  StDev    2.5%  97.5%      P
#> (Intercept)  0.2812 1.0151 -1.6637 2.3313 0.8010
#> sexfemale   -0.9065 0.5235 -1.9406 0.1066 0.0774
#> age          0.0136 0.0163 -0.0184 0.0453 0.3986
#> year         0.2462 0.0753  0.1030 0.3931 0.0004
#> 
#> Longitudinal outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev     2.5%   97.5% P
#> (Intercept) -8.2230 0.9798 -10.2546 -6.4654 0
#> year         0.4641 0.0683   0.3204  0.5958 0
#> age          0.0753 0.0164   0.0453  0.1096 0
#> 
#> MCMC summary:
#> chains: 1 
#> iterations per chain: 11000 
#> burn-in per chain: 1000 
#> thinning: 1 
#> time: 1.3 min

# }
```

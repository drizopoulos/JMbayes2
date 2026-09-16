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
#> marginal    4199.714 4985.926 -3083.326
#> conditional 3339.719 3172.861 -1812.772
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.9748 (Intr)
#> year   0.1767 0.3361
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%      P
#> drugD-penicil        -0.0255 0.2287 -0.4611 0.4265 0.9154
#> age                   0.0636 0.0091  0.0463 0.0818 0.0000
#> value(log(serBilir))  1.4162 0.1070  1.2227 1.6338 0.0000
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P
#> (Intercept)     0.2445 0.3564 -0.4556  0.9403 0.4928
#> year            0.2261 0.0363  0.1556  0.2973 0.0000
#> sexfemale      -0.2440 0.1804 -0.5965  0.1090 0.1746
#> I(year^2)       0.0026 0.0010  0.0007  0.0046 0.0002
#> age            -0.0017 0.0054 -0.0121  0.0090 0.7492
#> prothrombin     0.0529 0.0086  0.0366  0.0700 0.0000
#> year:sexfemale -0.0870 0.0380 -0.1615 -0.0135 0.0200
#> sigma           0.3457 0.0068  0.3323  0.3593 0.0000
#> 
#> MCMC summary:
#> chains: 1 
#> iterations per chain: 11000 
#> burn-in per chain: 1000 
#> thinning: 1 
#> time: 22 sec
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
#>                  DIC    WAIC      LPML
#> marginal    6635.479 6872.14 -3785.564
#> conditional 9087.043 8809.28 -4865.980
#> 
#> Random-effects covariance matrix:
#>                                                   
#>        StdDev   Corr                              
#> (Intr) 0.9922 (Intr)   year (Intr)    year (Intr) 
#> year   0.1743 0.3942                              
#> (Intr) 3.2229 0.5335 0.3416                       
#> year   0.5718 0.0592 0.3716 -0.3173               
#> (Intr) 2.8087 0.6176 0.4972 0.5285  0.0001        
#> year   0.4341 0.3977 0.6224 0.3794  0.2771 -0.0393
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%  97.5%      P
#> drugD-penicil        -0.1922 0.2621 -0.7099 0.3155 0.4564
#> age                   0.0326 0.0141  0.0020 0.0593 0.0350
#> value(log(serBilir))  0.6789 0.2203  0.2160 1.0936 0.0030
#> value(hepatomegaly)  -0.0412 0.0879 -0.2267 0.1258 0.6330
#> value(ascites)        0.5579 0.2077  0.1843 1.0865 0.0000
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P
#> (Intercept)     0.6729 0.1418  0.3950  0.9489 0.0000
#> year            0.2460 0.0281  0.1918  0.3017 0.0000
#> sexfemale      -0.2036 0.1448 -0.4901  0.0835 0.1570
#> year:sexfemale -0.0650 0.0295 -0.1236 -0.0087 0.0252
#> sigma           0.3481 0.0068  0.3354  0.3614 0.0000
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>                Mean  StDev    2.5%  97.5%      P
#> (Intercept)  0.0694 1.0134 -1.9184 2.0920 0.9532
#> sexfemale   -0.7647 0.5251 -1.8352 0.2462 0.1444
#> age          0.0147 0.0162 -0.0169 0.0461 0.3608
#> year         0.2541 0.0698  0.1190 0.3887 0.0002
#> 
#> Longitudinal outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev     2.5%   97.5% P
#> (Intercept) -8.9621 0.9847 -10.9888 -7.1858 0
#> year         0.5734 0.0681   0.4420  0.7205 0
#> age          0.0812 0.0160   0.0514  0.1147 0
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
#>                  DIC     WAIC      LPML
#> marginal    6637.878 7007.934 -3866.512
#> conditional 9035.769 8798.276 -4844.879
#> 
#> Random-effects covariance matrix:
#>                                                   
#>        StdDev   Corr                              
#> (Intr) 0.9913 (Intr)   year (Intr)    year  (Intr)
#> year   0.1804 0.4247                              
#> (Intr) 3.4220 0.5342 0.3458                       
#> year   0.6024 0.0343 0.3591 -0.3794               
#> (Intr) 2.3184 0.6455 0.5774 0.5968  -0.0947       
#> year   0.4166 0.4976 0.6947 0.4381  0.3218  0.3199
#> 
#> Survival outcome:
#>                         Mean  StDev    2.5%   97.5%      P
#> drugD-penicil        -0.1310 0.2646 -0.6481  0.3912 0.6074
#> age                   0.0476 0.0136  0.0198  0.0727 0.0074
#> value(log(serBilir))  0.9100 0.2021  0.5211  1.3053 0.0000
#> slope(log(serBilir))  3.5037 1.2792  1.2429  6.2379 0.0010
#> area(hepatomegaly)    0.1084 0.0811 -0.0502  0.2751 0.1652
#> value(ascites)       -0.5290 0.2552 -1.0677 -0.0867 0.0300
#> area(ascites)         0.8628 0.3374  0.3289  1.6461 0.0000
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P
#> (Intercept)     0.7042 0.1436  0.4273  0.9827 0.0000
#> year            0.2582 0.0297  0.2021  0.3187 0.0000
#> sexfemale      -0.2419 0.1478 -0.5352  0.0452 0.0996
#> year:sexfemale -0.0705 0.0295 -0.1312 -0.0149 0.0142
#> sigma           0.3481 0.0067  0.3351  0.3617 0.0000
#> 
#> Longitudinal outcome: hepatomegaly (family = binomial, link = logit)
#>                Mean  StDev    2.5%  97.5%      P
#> (Intercept)  0.2845 1.0364 -1.7434 2.3376 0.7930
#> sexfemale   -0.9183 0.5464 -1.9911 0.1413 0.0918
#> age          0.0138 0.0165 -0.0186 0.0461 0.4038
#> year         0.2545 0.0733  0.1101 0.3976 0.0000
#> 
#> Longitudinal outcome: ascites (family = binomial, link = logit)
#>                Mean  StDev     2.5%   97.5% P
#> (Intercept) -8.3055 0.9338 -10.2158 -6.6444 0
#> year         0.4722 0.0725   0.3256  0.6064 0
#> age          0.0761 0.0159   0.0471  0.1082 0
#> 
#> MCMC summary:
#> chains: 1 
#> iterations per chain: 11000 
#> burn-in per chain: 1000 
#> thinning: 1 
#> time: 1.2 min

# }
```

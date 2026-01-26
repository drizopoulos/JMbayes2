# Various Methods for Standard Generics

Methods for object of class `"jm"` for standard generic functions.

## Usage

``` r
coef(object, ...)

# S3 method for class 'jm'
coef(object, ...)

fixef(object, ...)

# S3 method for class 'jm'
fixef(object, outcome = Inf, ...)

ranef(object, ...)

# S3 method for class 'jm'
ranef(object, outcome = Inf, post_vars = FALSE, ...)

terms(x, ...)

# S3 method for class 'jm'
terms(x, process = c("longitudinal", "event"),
                      type = c("fixed", "random"), ...)

model.frame(formula, ...)

# S3 method for class 'jm'
model.frame(formula, process = c("longitudinal", "event"),
                            type = c("fixed", "random"), ...)

model.matrix(object, ...)

# S3 method for class 'jm'
model.matrix(object, ...)

family(object, ...)

# S3 method for class 'jm'
family(object, ...)

compare_jm(..., type = c("marginal", "conditional"),
  order = c("WAIC", "DIC", "LPML", "none"))
```

## Arguments

- object, x, formula:

  object inheriting from class `"jm"`.

- outcome:

  the index of the linear mixed submodel to extract the estimated fixed
  effects. If greater than the total number of submodels, extracts from
  all of them.

- post_vars:

  logical; if `TRUE`, returns the variance of the posterior
  distribution.

- process:

  which submodel(s) to extract the terms:  

  - if `"longitudinal"`, the linear mixed model(s), or

  - if `"event"`, the survival model.

- type:

  in `terms()` and `model.frame()`, which effects to select in the
  longitudinal process:  

  - if `"fixed"`, the fixed-effects, or

  - if `"random"`, the random-efects.

  in `compare_jm()`, which log-likelihood function use to calculate the
  criteria:  

  - if `"marginal"`, the marginal log-likelihood, or

  - if `"conditional"`, the conditional log-likelihood.

- ...:

  further arguments; currently, none is used.  
  in `compare_jm()`, a series of `jm` objects.

- order:

  which criteria use to sort the models in the output.

## Details

- `coef()`:

  Extracts estimated fixed effects for the event process from a fitted
  joint model.

- `fixef()`:

  Extracts estimated fixed effects for the longitudinal processes from a
  fitted joint model.

- `ranef()`:

  Extracts estimated random effects from a fitted joint model.

- `terms()`:

  Extracts the terms object(s) from a fitted joint model.

- `model.frame()`:

  Creates the model frame from a fitted joint model.

- `model.matrix()`:

  Creates the design matrices for linear mixed submodels from a fitted
  joint model.

- `family()`:

  Extracts the error distribution and link function used in the linear
  mixed submodel(s) from a fitted joint model.

- `compare_jm()`:

  Compares two or more fitted joint models using the criteria WAIC, DIC,
  and LPML.

## Value

- `coef()`:

  a list with the elements:  

  - `gammas`: estimated baseline fixed effects, and

  - `association`: estimated association parameters.

- `fixef()`:

  a numeric vector of the estimated fixed effects for the `outcome`
  selected. If the `outcome` is greater than the number of linear mixed
  submodels, it returns a list of numeric vectors for all outcomes.

- `ranef()`:

  a numeric matrix with rows denoting the individuals and columns the
  random effects. If `postVar = TRUE`, the numeric matrix has the extra
  attribute "postVar".

- `terms()`:

  if `process = "longitudinal"`, a list of the terms object(s) for the
  linear mixed model(s).  
  if `process = "event"`, the terms object for the survival model.

- `model.frame()`:

  if `process = "longitudinal"`, a list of the model frames used in the
  linear mixed model(s).  
  if `process = "event"`, the model frame used in the survival model.

- `model.matrix()`:

  a list of the design matrix(ces) for the linear mixed submodel(s).

- `family()`:

  a list of `family` objects.

- `compare_jm()`:

  a list with the elements:  

  - `table`: a table with the criteria calculated for each joint model,
    and

  - `type`: the log-likelihood function used to calculate the criteria.

## Author

Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>

## See also

[`jm`](https://drizopoulos.github.io/JMbayes2/reference/jm.md)

## Examples

``` r
# \donttest{
# linear mixed model fits
fit_lme1 <- lme(log(serBilir) ~ year:sex + age,
                random = ~ year | id, data = pbc2)

fit_lme2 <- lme(prothrombin ~ sex,
                random = ~ year | id, data = pbc2)

# cox model fit
fit_cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

# joint model fit
fit_jm <- jm(fit_cox, list(fit_lme1, fit_lme2), time_var = "year",
    n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)

# coef(): fixed effects for the event process
coef(fit_jm)
#> $gammas
#>       Mean 
#> 0.06070451 
#> 
#> $association
#> value(log(serBilir))   value(prothrombin) 
#>            1.3604613            0.1063424 
#> 

# fixef(): fixed effects for the first linear mixed submodel
fixef(fit_jm, outcome = 1)
#>    (Intercept)            age   year:sexmale year:sexfemale 
#>     0.66594783    -0.00445036     0.23584373     0.16253524 

# ranef(): random effects from all linear mixed submodels
head(ranef(fit_jm))
#>             [,1]        [,2]        [,3]         [,4]
#> [1,]  2.23319313  0.20276852  1.05194350  0.126098411
#> [2,] -0.36290322  0.00382308 -0.03130737  0.089480915
#> [3,] -0.19671871  0.06974686  0.51040201  0.210441867
#> [4,]  0.03344513  0.10239115  0.85942217  0.574548016
#> [5,]  0.32900876  0.22520488 -0.04018424  0.442201332
#> [6,] -0.62705384 -0.16161656 -0.11058573 -0.003544499

# terms(): random effects terms for the first linear mixed submodel
terms(fit_jm, process = "longitudinal", type = "random")[[1]]
#> ~year
#> attr(,"variables")
#> list(year)
#> attr(,"factors")
#>      year
#> year    1
#> attr(,"term.labels")
#> [1] "year"
#> attr(,"order")
#> [1] 1
#> attr(,"intercept")
#> [1] 1
#> attr(,"response")
#> [1] 0
#> attr(,".Environment")
#> <environment: R_GlobalEnv>
#> attr(,"predvars")
#> list(year)
#> attr(,"dataClasses")
#>      year 
#> "numeric" 

# mode.frame(): model frame for the fixed effects in the second
# linear mixed submodel
head(model.frame(fit_jm, process = "longitudinal", type = "fixed")[[2]])
#>   prothrombin    sex
#> 1        12.2 female
#> 2        11.2 female
#> 3        10.6 female
#> 4        11.0 female
#> 5        11.6 female
#> 6        10.6 female

# model.matrix(): fixed effects design matrix for the first linear
# mixed submodel
head(model.matrix(fit_jm)[[1]])
#>   (Intercept)      age year:sexmale year:sexfemale
#> 1           1 58.76684            0      0.0000000
#> 2           1 58.76684            0      0.5256817
#> 3           1 56.44782            0      0.0000000
#> 4           1 56.44782            0      0.4983025
#> 5           1 56.44782            0      0.9993429
#> 6           1 56.44782            0      2.1027270

# family(): family objects from both linear mixed submodels
family(fit_jm)
#> [[1]]
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> 
#> [[2]]
#> 
#> Family: gaussian 
#> Link function: identity 
#> 
#> 

# compare_jm(): compare two fitted joint models
fit_lme1b <- lme(log(serBilir) ~ 1,
                  random = ~ year | id, data = pbc2)

fit_jm2 <- jm(fit_cox, list(fit_lme1b, fit_lme2), time_var = "year",
    n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)

compare_jm(fit_jm, fit_jm2)
#> 
#>               DIC     WAIC      LPML
#>  fit_jm2 10512.30 10540.46 -5268.228
#>   fit_jm 10665.69 11146.24 -6097.345
#> 
#> The criteria are calculated based on the marginal log-likelihood.
# }
```

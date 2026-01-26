# Slice-aware model fitting generics

These functions are S3 generics exported by JMbayes2 to enable method
dispatch on the class of the `data` argument. When `data` is a
`"sliced_data"` object (as produced by
[`slicer`](https://drizopoulos.github.io/JMbayes2/reference/slicer.md)),
the corresponding slice-wise fitting methods are used.

## Usage

``` r
lme(fixed, data, ...)
coxph(formula, data, ...)
mixed_model(fixed, data, ...)
```

## Arguments

- fixed:

  For `lme()` and `mixed_model()`, a model formula for the fixed
  effects.

- formula:

  For `coxph()`, a survival model formula.

- data:

  A data frame (default methods) or a `"sliced_data"` object (slice-wise
  methods).

- ...:

  Further arguments passed to the underlying model-fitting functions.
  For slice-wise methods, additional arguments include `parallel_out`
  and `cores` to control parallel execution across slices.

## Details

When `data` is a regular data frame, the default methods call
[`nlme::lme()`](https://rdrr.io/pkg/nlme/man/lme.html),
[`survival::coxph()`](https://rdrr.io/pkg/survival/man/coxph.html), and
[`GLMMadaptive::mixed_model()`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.html),
respectively.

When `data` is a `"sliced_data"` object, the corresponding
`*.sliced_data` methods fit the requested model independently within
each slice and return a list of fits (one per slice).

## Note

JMbayes2 exports S3 generics `lme()`, `coxph()`, and `mixed_model()` to
enable dispatch on `"sliced_data"`. When JMbayes2 is attached, these
names mask [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html),
[`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html), and
[`GLMMadaptive::mixed_model`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.html).
Use the `pkg::fun` form to call the original functions.

## See also

[`slicer`](https://drizopoulos.github.io/JMbayes2/reference/slicer.md),
[`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html),
[`survival::coxph`](https://rdrr.io/pkg/survival/man/coxph.html),
[`GLMMadaptive::mixed_model`](https://drizopoulos.github.io/GLMMadaptive/reference/mixed_model.html).

## Examples

``` r
if (FALSE) { # \dontrun{

slc <- slicer(n_slices = 2, id_var = "id", data_long = pbc2, data_surv = pbc2.id)
n_cores <- max(parallel::detectCores() - 1L, 1L)

lme_fit <- lme(fixed  = log(serBilir) ~ year * sex,
               data   = slc$long,
               random = ~ year | id,
               cores  = n_cores)

cox_fit <- coxph(formula = Surv(years, status2) ~ sex,
                 data    = slc$surv,
                 cores   = n_cores)

mxm_fit <- mixed_model(fixed  = ascites ~ year + sex,
                       data   = slc$long,
                       random = ~ year | id,
                       family = binomial(),
                       cores  = n_cores)
} # }
```

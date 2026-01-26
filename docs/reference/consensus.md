# Consensus Combination of Posterior Draws from Joint Model Fits

This function combines posterior draws from multiple joint model fits
into a single set of draws.

## Usage

``` r
consensus(object, parm, method = c("union", "equal_weight", "var_weight"),
          seed = 123L)
```

## Arguments

- object:

  an object of class `"sliced_jm"`.

- parm:

  a character vector with the names of parameters to combine. These must
  correspond to elements in `fit$mcmc[[parm]]` for each fit (e.g.,
  `"gammas"`, `"alphas"`, `"betas1"`, `"betas2"`).

- method:

  the consensus method used to combine draws.

  `"union"`

  :   concatenate draws across slices (no averaging).

  `"equal_weight"`

  :   compute an iteration-wise simple average across slices.

  `"var_weight"`

  :   compute an iteration-wise weighted average across slices using
      inverse-variance weights.

- seed:

  an integer seed used for the within-slice random permutation step
  (used by `"equal_weight"` and `"var_weight"`).

## Value

An object of class `"consensus_jm"` with components:

- `method`:

  the selected method.

- `parm`:

  the requested parameter(s) block(s).

- `n_splits`:

  number of fits combined.

- `draws`:

  named list of combined draw matrices (one per parameter block).

- `weights`:

  for `"equal_weight"` and `"var_weight"`, a named list of weight
  matrices.

- `n_draws`:

  number of combined draws per parameter block.

- `summary`:

  named list of summary matrices with `Mean`, `StDev`, `2.5%`, `97.5%`,
  and `P`).

- `seed`:

  seed used.

## Author

Pedro Miranda-Afonso <p.mirandaafonso@erasmusmc.nl>

# Variogram for Longitudinal Data

It computes the semi-variogram for longitudinal data.

## Usage

``` r
variogram(y, times, id)
```

## Arguments

- y:

  a numeric vector of longitudinal responses.

- times:

  a numeric vector of times at which the longitudinal responses were
  collected.

- id:

  a numeric vector of a factor of subject id numbers.

## Value

A list with two components, i.e., `svar`, a two-column matrix with the
time lags and the variogram values, and `sigma2` the total variance.

## Author

Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>

## Examples

``` r
ind <- aids$patient == 2
yy <- aids$CD4[ind]
tt <- aids$obstime[ind]
ids <- aids$patient[ind]
variogram(yy, tt, ids)
#> $svar
#>      time_lag     diffs2
#> [1,]        6 1.61906969
#> [2,]       12 1.51724651
#> [3,]       18 0.87722340
#> [4,]        6 6.27097906
#> [5,]       12 4.87980798
#> [6,]        6 0.08712153
#> 
#> $sigma2
#> [1] NA
#> 
#> attr(,"class")
#> [1] "vrgm"
```

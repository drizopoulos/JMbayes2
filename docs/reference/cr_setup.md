# Transform Competing Risks Data in Long Format

In a competing risks setting this function expands the data frame with a
single row per subject to a data frame in the long format in which each
subject has as many rows as the number of competing events.

## Usage

``` r
crisk_setup(data, statusVar, censLevel,
    nameStrata = "strata", nameStatus = "status2")
```

## Arguments

- data:

  the data frame containing the competing risk data with a single row
  per subject.

- statusVar:

  a character string denoting the name of the variable in `data` that
  identifies the status variable which equals 1 if the subject had any
  of the competing events and 0 otherwise.

- censLevel:

  a character string or a scalar denoting the censoring level in the
  `statusVar` variable of `data`.

- nameStrata:

  a character string denoting the variable that will be added in the
  long version of `data` denoting the various causes of event.

- nameStatus:

  a character string denoting the variable that will be added in the
  long version of `data` denoting if the subject experience any of the
  competing events.

## Value

A data frame in the long format with multiple rows per subject.

## References

Rizopoulos, D. (2012) *Joint Models for Longitudinal and Time-to-Event
Data: with Applications in R*. Boca Raton: Chapman and Hall/CRC.

Putter, H., Fiocco, M., and Geskus, R. (2007). Tutorial in
biostatistics: Competing risks and multi-state models. *Statistics in
Medicine* **26**, 2389–2430.

## Author

Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>

## Examples

``` r
head(crisk_setup(pbc2.id, "status", "alive"))
#>     id     years status      drug      age    sex year ascites hepatomegaly
#> 1    1  1.095170   dead D-penicil 58.76684 female    0     Yes          Yes
#> 1.1  1  1.095170   dead D-penicil 58.76684 female    0     Yes          Yes
#> 2    2 14.152338  alive D-penicil 56.44782 female    0      No          Yes
#> 2.1  2 14.152338  alive D-penicil 56.44782 female    0      No          Yes
#> 3    3  2.770781   dead D-penicil 70.07447   male    0      No           No
#> 3.1  3  2.770781   dead D-penicil 70.07447   male    0      No           No
#>     spiders                   edema serBilir serChol albumin alkaline  SGOT
#> 1       Yes edema despite diuretics     14.5     261    2.60     1718 138.0
#> 1.1     Yes edema despite diuretics     14.5     261    2.60     1718 138.0
#> 2       Yes                No edema      1.1     302    4.14     7395 113.5
#> 2.1     Yes                No edema      1.1     302    4.14     7395 113.5
#> 3        No      edema no diuretics      1.4     176    3.48      516  96.1
#> 3.1      No      edema no diuretics      1.4     176    3.48      516  96.1
#>     platelets prothrombin histologic status2       strata
#> 1         190        12.2          4       1         dead
#> 1.1       190        12.2          4       0 transplanted
#> 2         221        10.6          3       0         dead
#> 2.1       221        10.6          3       0 transplanted
#> 3         151        12.0          4       1         dead
#> 3.1       151        12.0          4       0 transplanted
```

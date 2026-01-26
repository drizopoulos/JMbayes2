# Didanosine versus Zalcitabine in HIV Patients

A randomized clinical trial in which both longitudinal and survival data
were collected to compare the efficacy and safety of two antiretroviral
drugs in treating patients who had failed or were intolerant of
zidovudine (AZT) therapy.

## Format

A data frame with 1408 observations on the following 9 variables.

- `patient`:

  patients identifier; in total there are 467 patients.

- `Time`:

  the time to death or censoring.

- `death`:

  a numeric vector with 0 denoting censoring and 1 death.

- `CD4`:

  the CD4 cells count.

- `obstime`:

  the time points at which the CD4 cells count was recorded.

- `drug`:

  a factor with levels `ddC` denoting zalcitabine and `ddI` denoting
  didanosine.

- `gender`:

  a factor with levels `female` and `male`.

- `prevOI`:

  a factor with levels `AIDS` denoting previous opportunistic infection
  (AIDS diagnosis) at study entry, and `noAIDS` denoting no previous
  infection.

- `AZT`:

  a factor with levels `intolerance` and `failure` denoting AZT
  intolerance and AZT failure, respectively.

## References

Goldman, A., Carlin, B., Crane, L., Launer, C., Korvick, J., Deyton, L.
and Abrams, D. (1996) Response of CD4+ and clinical consequences to
treatment using ddI or ddC in patients with advanced HIV infection.
*Journal of Acquired Immune Deficiency Syndromes and Human
Retrovirology* **11**, 161–169.

Guo, X. and Carlin, B. (2004) Separate and joint modeling of
longitudinal and event time data using standard computer packages. *The
American Statistician* **58**, 16–24.

## Note

The data frame `aids.id` contains the first CD4 cell count measurement
for each patient. This data frame is used to fit the survival model.

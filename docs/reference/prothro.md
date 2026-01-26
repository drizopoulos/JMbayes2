# Prednisone versus Placebo in Liver Cirrhosis Patients

A randomized trial on 488 liver cirrhosis patients.

## Format

Two data frames with the following variables.

- `id`:

  patients identifier; in total there are 467 patients.

- `pro`:

  prothrobin measurements.

- `time`:

  for data frame `prothro` the time points at which the prothrobin
  measurements were taken; for data frame `prothros` the time to death
  or censoring.

- `death`:

  a numeric vector with 0 denoting censoring and 1 death.

- `treat`:

  randomized treatment; a factor with levels "placebo" and "prednisone".

## Source

<http://www.gllamm.org/books/readme.html#14.6>.

## References

Andersen, P. K., Borgan, O., Gill, R. D. and Keiding, N. (1993).
*Statistical Models Based on Counting Processes*. New York: Springer.

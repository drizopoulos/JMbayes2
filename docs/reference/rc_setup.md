# Combine Recurring and Terminal Event Data in Long Format

This function combines two data frames, the recurring-event and
terminal-event/competing-risks datasets, into one. Each subject has as
many rows in the new data frame as the number of recurrent risk periods
plus one for each terminal event/competing risk.

## Usage

``` r
rc_setup(rc_data, trm_data,
    idVar = "id", statusVar = "status",
    startVar = "start", stopVar = "stop",
    trm_censLevel,
    nameStrata = "strata", nameStatus = "status")
```

## Arguments

- rc_data:

  the data frame containing the recurring-event data with multiple rows
  per subject.

- trm_data:

  the data frame containing the terminal-event/competing-risks data with
  a single row per subject.

- idVar:

  a character scalar denoting the name of the variable in `rc_data` and
  `trm_data` that identifies the subject/group.

- statusVar:

  a character string denoting the name of the variable in `rc_data` and
  `trm_data` that identifies the status variable. In `rc_data` equals 1
  if the subject had an event and 0 otherwise. In `trm_data` equals to
  the event or censoring level.

- startVar:

  a character string denoting the name of the variable in `rc_data` that
  identifies the starting time for the risk interval.

- stopVar:

  a character string denoting the name of the variable in `rc_data` and
  `trm_data` that identifies the event or censoring time.

- trm_censLevel:

  a character string or a scalar denoting the censoring level in the
  statusVar variable of `trm_data`.

- nameStrata:

  a character string denoting the variable that will be added in the
  long version of `data` denoting the various causes of event.

- nameStatus:

  a character string denoting the variable that will be added in the
  long version of `data` denoting if the subject had an event.

## Value

A data frame in the long format with multiple rows per subject.

## Author

Pedro Miranda-Afonso <p.mirandaafonso@erasmusmc.nl>

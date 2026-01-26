# Mayo Clinic Primary Biliary Cirrhosis Data

Follow up of 312 randomised patients with primary biliary cirrhosis, a
rare autoimmune liver disease, at Mayo Clinic.

## Format

A data frame with 1945 observations on the following 20 variables.

- `id`:

  patients identifier; in total there are 312 patients.

- `years`:

  number of years between registration and the earlier of death,
  transplantion, or study analysis time.

- `status`:

  a factor with levels `alive`, `transplanted` and `dead`.

- `drug`:

  a factor with levels `placebo` and `D-penicil`.

- `age`:

  at registration in years.

- `sex`:

  a factor with levels `male` and `female`.

- `year`:

  number of years between enrollment and this visit date, remaining
  values on the line of data refer to this visit.

- `ascites`:

  a factor with levels `No` and `Yes`.

- `hepatomegaly`:

  a factor with levels `No` and `Yes`.

- `spiders`:

  a factor with levels `No` and `Yes`.

- `edema`:

  a factor with levels `No edema` (i.e., no edema and no diuretic
  therapy for edema), `edema no diuretics` (i.e., edema present without
  diuretics, or edema resolved by diuretics), and
  `edema despite diuretics` (i.e., edema despite diuretic therapy).

- `serBilir`:

  serum bilirubin in mg/dl.

- `serChol`:

  serum cholesterol in mg/dl.

- `albumin`:

  albumin in g/dl.

- `alkaline`:

  alkaline phosphatase in U/liter.

- `SGOT`:

  SGOT in U/ml.

- `platelets`:

  platelets per cubic ml / 1000.

- `prothrombin`:

  prothrombin time in seconds.

- `histologic`:

  histologic stage of disease.

- `status2`:

  a numeric vector with the value 1 denoting if the patient was dead,
  and 0 if the patient was alive or transplanted.

## References

Fleming, T. and Harrington, D. (1991) *Counting Processes and Survival
Analysis*. Wiley, New York.

Therneau, T. and Grambsch, P. (2000) *Modeling Survival Data: Extending
the Cox Model*. Springer-Verlag, New York.

## Note

The data frame `pbc2.id` contains the first measurement for each
patient. This data frame is used to fit the survival model.

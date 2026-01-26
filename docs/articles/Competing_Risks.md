# Competing Risks

## Competing Risks

### Prepare data

The first step in fitting a joint model for competing events in
**JMbayes2** is to prepare the data for the event process. If there are
K competing events, each subject must have K rows, one for each possible
cause. The observed event time T_i of each subject is repeated K times,
and there are two indicator variables, namely one identifying the cause
and one indicating whether the corresponding event type is the one that
occurred. Standard survival datasets that include a single row per
patient can be easily transformed to the competing risks long format
using the function
[`crisk_setup()`](https://drizopoulos.github.io/JMbayes2/reference/cr_setup.md).
This function accepts as main arguments the survival data in the
standard format with a single row per patient, the name of the status
variable, and the level in this status variable that corresponds to
censoring. We illustrate the use of this function in the PBC data, where
we treat as competing risks transplantation and death:

``` r
pbc2.id[pbc2.id$id %in% c(1, 2, 5), c("id", "years", "status")]
#>   id     years       status
#> 1  1  1.095170         dead
#> 2  2 14.152338        alive
#> 5  5  4.120578 transplanted

pbc2.idCR <- crisk_setup(pbc2.id, statusVar = "status", censLevel = "alive", 
                         nameStrata = "CR")

pbc2.idCR[pbc2.idCR$id %in% c(1, 2, 5), 
          c("id", "years", "status", "status2", "CR")]
#>     id     years       status status2           CR
#> 1    1  1.095170         dead       1         dead
#> 1.1  1  1.095170         dead       0 transplanted
#> 2    2 14.152338        alive       0         dead
#> 2.1  2 14.152338        alive       0 transplanted
#> 5    5  4.120578 transplanted       0         dead
#> 5.1  5  4.120578 transplanted       1 transplanted
```

Note that each patient is now represented by two rows (we have two
possible causes of discontinuation from the study, death, and
transplantation), the event time variable `years` is identical in both
rows of each patient, variable `CR` denotes the cause for the specific
line of the long dataset, and variable `status2` equals 1 if the
corresponding event occurred.

### Fit models

For the event process, we specify cause-specific relative risk models.
Using dataset `pbc2.idCR`, we fit the corresponding cause-specific Cox
regressions by including the interaction terms of age and treatment with
variable `CR`, which is treated as a stratification variable using the
`strata()` function:

``` r
CoxFit_CR <- coxph(Surv(years, status2) ~ (age + drug):strata(CR),
                     data = pbc2.idCR)
```

We include two longitudinal outcomes for the longitudinal process: serum
bilirubin and the prothrombin time. For the former, we use quadratic
orthogonal polynomials in the fixed- and random-effects parts, and for
the latter, linear evolutions:

``` r
fm1 <- lme(log(serBilir) ~ poly(year, 2) * drug, data = pbc2, 
           random = ~ poly(year, 2) | id)
fm2 <- lme(prothrombin ~ year * drug, data = pbc2, random = ~ year | id)
```

To specify that each longitudinal outcome has a separate association
coefficient per competing risk, we define the corresponding functional
forms:

``` r
CR_forms <- list(
    "log(serBilir)" = ~ value(log(serBilir)):CR,
    "prothrombin" = ~ value(prothrombin):CR
)
```

Finally, the competing risks joint model is fitted with the following
call to [`jm()`](https://drizopoulos.github.io/JMbayes2/reference/jm.md)
(due to the complexity of the model, we have increased the number of
MCMC iterations and the burn-in period per chain). Also, because
relatively few patients received a transplantation, we specify a Weibull
baseline hazard function for this competing event, and the default
penalized B-spline approximation for death:

``` r
jFit_CR <- jm(CoxFit_CR, list(fm1, fm2), time_var = "year", 
              functional_forms = CR_forms, 
              baseline_hazard = c("weibull", NA),
              n_iter = 25000L, n_burnin = 5000L, n_thin = 5L)
#> Warning in jm.default(CoxFit_CR, list(fm1, fm2), time_var = "year",
#> functional_forms = CR_forms, : unknown names in control: baseline_hazard

summary(jFit_CR)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = CoxFit_CR, Mixed_objects = list(fm1, 
#>     fm2), time_var = "year", functional_forms = CR_forms, baseline_hazard = c("weibull", 
#>     NA), n_iter = 25000L, n_burnin = 5000L, n_thin = 5L)
#> 
#> Data Descriptives:
#> Number of Groups: 312        Number of events: 169 (27.1%)
#> Number of Observations:
#>   log(serBilir): 1945
#>   prothrombin: 1945
#> 
#>                  DIC     WAIC      LPML
#> marginal    10824.49 11528.31 -6753.594
#> conditional 15751.80 15442.25 -8231.393
#> 
#> Random-effects covariance matrix:
#>                                              
#>        StdDev    Corr                        
#> (Intr) 1.3404  (Intr)  p(,2)1  p(,2)2  (Intr)
#> p(,2)1 23.0867 0.6987                        
#> p(,2)2 12.4129 -0.2690 -0.1464               
#> (Intr) 0.7873  0.6343  0.4354  -0.3299       
#> year   0.3269  0.4310  0.3351  -0.0589 0.0307
#> 
#> Survival Outcome:
#>                                         Mean  StDev    2.5%   97.5%      P
#> age:strata(CR)transplanted           -0.0711 0.0299 -0.1299 -0.0137 0.0130
#> age:strata(CR)dead                    0.0632 0.0100  0.0435  0.0832 0.0000
#> drugD-penicil:strata(CR)transplanted -0.2895 0.4166 -1.1253  0.5072 0.4852
#> drugD-penicil:strata(CR)dead         -0.0067 0.1873 -0.3681  0.3639 0.9710
#> value(log(serBilir)):CRtransplanted   1.0296 0.2206  0.6229  1.4844 0.0000
#> value(log(serBilir)):CRdead           1.4560 0.1169  1.2343  1.6939 0.0000
#> value(prothrombin):CRtransplanted     0.1262 0.1519 -0.1746  0.4012 0.4428
#> value(prothrombin):CRdead             0.1387 0.0509  0.0280  0.2304 0.0125
#>                                        Rhat
#> age:strata(CR)transplanted           1.0261
#> age:strata(CR)dead                   1.0144
#> drugD-penicil:strata(CR)transplanted 1.0149
#> drugD-penicil:strata(CR)dead         1.0015
#> value(log(serBilir)):CRtransplanted  1.0172
#> value(log(serBilir)):CRdead          1.0224
#> value(prothrombin):CRtransplanted    1.3322
#> value(prothrombin):CRdead            1.0244
#> 
#> Longitudinal Outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev     2.5%   97.5%      P   Rhat
#> (Intercept)     1.1949 0.1139   0.9710  1.4193 0.0000 0.9998
#> poly(year, 2)1 27.7281 3.0802  21.9088 34.1392 0.0000 1.0096
#> poly(year, 2)2  1.1796 1.7760  -2.3127  4.6065 0.5030 1.0031
#> drugD-penicil  -0.1902 0.1578  -0.5039  0.1232 0.2262 1.0001
#> p(,2)1         -3.2070 3.6367 -10.4855  3.8662 0.3783 1.0031
#> p(,2)2         -1.0882 2.1796  -5.2804  3.1996 0.6177 1.0027
#> sigma           0.3023 0.0062   0.2903  0.3147 0.0000 1.0008
#> 
#> Longitudinal Outcome: prothrombin (family = gaussian, link = identity)
#>                       Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)        10.6331 0.0829 10.4704 10.7984 0.0000 1.0029
#> year                0.2944 0.0395  0.2178  0.3719 0.0000 1.0029
#> drugD-penicil      -0.0954 0.1165 -0.3215  0.1299 0.4148 1.0007
#> year:drugD-penicil -0.0238 0.0518 -0.1261  0.0780 0.6453 1.0005
#> sigma               1.0556 0.0204  1.0166  1.0971 0.0000 1.0013
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 25000 
#> burn-in per chain: 5000 
#> thinning: 5 
#> time: 7 min
```

### Dynamic predictions

Based on the fitted competing risks joint model, we will illustrate how
(dynamic) predictions can be calculated for the cause-specific
cumulative risk probabilities. As an example, we will show these
calculations for Patient 81 from the PBC dataset. First, we extract the
data on this subject.

``` r
ND_long <- pbc2[pbc2$id == 81, ]
ND_event <- pbc2.idCR[pbc2.idCR$id == 81, ]
ND_event$status2 <- 0
ND <- list(newdataL = ND_long, newdataE = ND_event)
```

The first line extracts the longitudinal measurements, and the second
line extracts the event times per cause (i.e., death and
transplantation). This patient died at 6.95 years, but to make the
calculation of cause-specific cumulative risk more relevant, we presume
that she did not have the event, and we set the event status variable
`status2` to zero. The last line combines the two datasets in a list.
*Note:* this last step is a prerequisite from the
[`predict()`](https://rdrr.io/r/stats/predict.html) method for competing
risks joint model. That is, the datasets provided in the arguments
`newdata` and `newdata2` need to be named lists with two components. The
first component needs to be named `newdataL` and contain the dataset
with the longitudinal measurements. The second component needs to be
named `newdataE` and contain the dataset with the event information.

The predictions are calculated using the
[`predict()`](https://rdrr.io/r/stats/predict.html) method. The first
call to this function calculates the prediction for the longitudinal
outcomes at the times provided in the `times` argument, and the second
call calculates the cause-specific cumulative risk probabilities. By
setting the argument `return_newdata` to `TRUE` in both calls, we can
use the corresponding
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) method to
depict the predictions:

``` r
predLong <- predict(jFit_CR, newdata = ND, return_newdata = TRUE,
                    times = seq(6.5, 15, length = 25))

predEvent <- predict(jFit_CR, newdata = ND, return_newdata = TRUE,
                     process = "event")

plot(predLong, predEvent, outcomes = 1:2, ylim_long_outcome_range = FALSE,
     col_line_event = c("#03BF3D", "#FF0000"), 
     fill_CI_event = c("#03BF3D4D", "#FF00004D"), pos_ylab_long = c(1.5, 11.5))
legend(x = 8.1, y = 0.45, legend = levels(pbc2.idCR$CR), 
       lty = 1, lwd = 2, col = c("#03BF3D", "#FF0000"), bty = "n", cex = 0.8)
```

![](Competing_Risks_files/figure-html/CIFs-1.png)

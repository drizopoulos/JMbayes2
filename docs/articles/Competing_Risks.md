# Competing Risks

## Joint Models with Competing Risks

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
              base_hazard = c("weibull", NA),
              n_iter = 25000L, n_burnin = 5000L, n_thin = 5L)

summary(jFit_CR)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = CoxFit_CR, Mixed_objects = list(fm1, 
#>     fm2), time_var = "year", functional_forms = CR_forms, base_hazard = c("weibull", 
#>     NA), n_iter = 25000L, n_burnin = 5000L, n_thin = 5L)
#> 
#> Data Descriptives:
#> Number of groups: 312        Number of events: 169 (27.1%)
#> Number of observations:
#>   log(serBilir): 1945
#>   prothrombin: 1945
#> 
#>                  DIC     WAIC      LPML
#> marginal    10836.48 11593.33 -6844.690
#> conditional 15758.81 15447.63 -8243.271
#> 
#> Random-effects covariance matrix:
#>                                              
#>        StdDev    Corr                        
#> (Intr) 1.3494  (Intr)  p(,2)1  p(,2)2  (Intr)
#> p(,2)1 23.2895 0.7128                        
#> p(,2)2 12.3887 -0.2599 -0.1483               
#> (Intr) 0.7892  0.6336  0.4420  -0.3325       
#> year   0.3275  0.4354  0.3442  -0.0509 0.0395
#> 
#> Survival outcome:
#>                                         Mean  StDev    2.5%   97.5%      P
#> age:strata(CR)transplanted           -0.0789 0.0256 -0.1303 -0.0314 0.0000
#> age:strata(CR)dead                    0.0649 0.0099  0.0456  0.0844 0.0000
#> drugD-penicil:strata(CR)transplanted -0.2747 0.3968 -1.0631  0.4999 0.4957
#> drugD-penicil:strata(CR)dead          0.0043 0.1882 -0.3631  0.3736 0.9913
#> value(log(serBilir)):CRtransplanted   1.0493 0.2150  0.6449  1.4885 0.0000
#> value(log(serBilir)):CRdead           1.4625 0.1164  1.2446  1.6966 0.0000
#> value(prothrombin):CRtransplanted    -0.0093 0.1484 -0.3118  0.2625 0.9657
#> value(prothrombin):CRdead             0.1520 0.0461  0.0593  0.2381 0.0025
#>                                        Rhat
#> age:strata(CR)transplanted           1.0439
#> age:strata(CR)dead                   1.0020
#> drugD-penicil:strata(CR)transplanted 1.0029
#> drugD-penicil:strata(CR)dead         1.0028
#> value(log(serBilir)):CRtransplanted  1.0039
#> value(log(serBilir)):CRdead          1.0033
#> value(prothrombin):CRtransplanted    1.0387
#> value(prothrombin):CRdead            1.0189
#> 
#> Longitudinal outcome: log(serBilir) (family = gaussian, link = identity)
#>                   Mean  StDev     2.5%   97.5%      P   Rhat
#> (Intercept)     1.2020 0.1154   0.9794  1.4270 0.0000 1.0012
#> poly(year, 2)1 28.0623 3.0988  22.2541 34.4407 0.0000 1.0047
#> poly(year, 2)2  1.2198 1.7764  -2.2226  4.6556 0.4922 1.0023
#> drugD-penicil  -0.1915 0.1583  -0.4986  0.1238 0.2240 1.0011
#> p(,2)1         -3.2039 3.6115 -10.3220  3.7644 0.3695 1.0096
#> p(,2)2         -1.0293 2.1549  -5.2776  3.2148 0.6235 1.0047
#> sigma           0.3024 0.0062   0.2906  0.3150 0.0000 1.0015
#> 
#> Longitudinal outcome: prothrombin (family = gaussian, link = identity)
#>                       Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)        10.6346 0.0834 10.4711 10.8008 0.0000 1.0013
#> year                0.2936 0.0395  0.2190  0.3725 0.0000 1.0024
#> drugD-penicil      -0.0948 0.1173 -0.3279  0.1376 0.4172 1.0002
#> year:drugD-penicil -0.0231 0.0516 -0.1255  0.0778 0.6558 1.0005
#> sigma               1.0547 0.0204  1.0153  1.0956 0.0000 1.0001
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 25000 
#> burn-in per chain: 5000 
#> thinning: 5 
#> time: 2.9 min
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

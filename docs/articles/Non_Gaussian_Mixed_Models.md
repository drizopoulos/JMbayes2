# Non-Gaussian Mixed Models

## Non-Gaussian Joint Models with JMbayes2

Taking advantage of the versatility of the
[**GLMMadaptive**](https://drizopoulos.github.io/GLMMadaptive/) package,
**JMbayes2** can fit joint models with several different types of
mixed-effects models. The following examples illustrate these
capabilities. All examples have the same structure, namely, first, a
short motivation for each mixed-model is given, followed by a piece of R
code simulating data from a joint model with the respective
mixed-effects sub-model, closing by the syntax to fit the joint model.
In this last part, the main difference per example is the call to
[`mixed_model()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md).

### Beta mixed models

With very few exceptions, continuous outcomes that we wish to analyze
have some natural bounds. For example, the levels of blood biomarkers
for a set of patients. However, most observations are often located far
away from these natural bounds, and an assumption of a normal
distribution for the outcome can be safely made. In some settings,
though, we can have outcomes for which a substantial percentage of the
observations are located near the boundaries, leading to skewed or
U-shaped distributions. A linear mixed model with normal error terms
often does not fit such longitudinal outcomes well. A natural
alternative is to select a distribution that respects the bounded nature
of the outcome. The most well-known distribution for such outcomes is
the Beta distribution defined in the (0, 1) interval (*note:* a bounded
outcome Y^\* in the (a, b) interval can be transformed to the Y =
(Y^\* - a) / (b - a) in the (0, 1) interval).

The following code illustrates how to simulate data from a joint model
with a Beta mixed effects model. The default functional form is assumed,
i.e., that the linear predictor \eta(t) of the mixed model is associated
with the hazard of an event at time t. The linear predictor is related
to the mean \mu(t) of the Beta distribution under the logit link
function, i.e., \log\[\mu(t) / \\1 - \mu(t)\\\] = \eta(t).

``` r

set.seed(1234)
n <- 200 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 10 # maximum follow-up time

# we construct a data frame with the design:
# everyone has a baseline measurement, and then measurements at random 
# follow-up times up to t_max
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ time, data = DF)

betas <- c(-2.2, -0.25, 0.24, -0.05) # fixed effects coefficients
phi <- 5 # precision parameter of the Beta distribution
D11 <- 1.0 # variance of random intercepts
D22 <- 0.5 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
# mean of the Beta distribution
mu_y <- plogis(eta_y) # plogis(eta_y) = exp(eta_y) / (1 + exp(eta_y))
# we simulate Beta longitudinal data
DF$y <- rbeta(n * K, shape1 = mu_y * phi, shape2 = phi * (1 - mu_y))
# we transform to (0, 1)
DF$y <- (DF$y * (nrow(DF) - 1) + 0.5) / nrow(DF)

upp_Cens <- 15 # fixed Type I censoring time
shape_wb <- 5 # shape Weibull
alpha <- 0.8 # association coefficients
gammas <- c("(Intercept)" = -9, "sex" = 0.5)
W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
# linear predictor for the survival model
eta_t <- as.vector(W %*% gammas)
# to simulate event times we use inverse transform sampling
# (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
# to find t, such that S(t) = u, where S(.) is the survival function, and u a 
# number from the Unif(0, 1) distribution. The function below calculates 
# log(u) - log(S(t)), and for a given u, we want to find t for which it equals
# zero. We do that below using the uniroot() function
invS <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h <- function (s) {
    X_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f <- as.vector(X_at_s %*% betas +
                     rowSums(Z_at_s * b[rep(i, nrow(Z_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f * alpha)
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h, lower = 0, upper = t)$value + log(u[i])
}
# we simulate the event times
u <- runif(n)
trueTimes <- numeric(n)
for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}

# we use fixed Type I right censoring denoting the end of the trial.
Ctimes <- upp_Cens
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator

# we keep the longitudinal measurements before the event times
DF$Time <- Time[DF$id]
DF$event <- event[DF$id]
DF <- DF[DF$time <= DF$Time, ]
```

To fit the corresponding joint model, we fit first a Beta mixed model
using the
[`beta.fam()`](https://rdrr.io/pkg/GLMMadaptive/man/extra_fams.html)
family object into the call of
[`mixed_model()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md):

``` r

DF_id <- DF[!duplicated(DF$id), ]
Cox_fit <- coxph(Surv(Time, event) ~ sex, data = DF_id)
Beta_MixMod <- mixed_model(y ~ sex * time, random = ~ time | id, data = DF,
                           family = beta.fam())

jointFit <- jm(Cox_fit, Beta_MixMod, time_var = "time")
summary(jointFit)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = Cox_fit, Mixed_objects = Beta_MixMod, 
#>     time_var = "time")
#> 
#> Data Descriptives:
#> Number of groups: 200        Number of events: 158 (79%)
#> Number of observations:
#>   y: 1182
#> 
#>                   DIC      WAIC     LPML
#> marginal    -4014.132 -4054.537 1871.975
#> conditional -3827.281 -3930.295 1806.311
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.8555 (Intr)
#> time   0.4712 0.0696
#> 
#> Survival outcome:
#>             Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale 0.2219 0.3010 -0.3670 0.8190 0.4647 1.0125
#> value(y)  1.0608 0.1103  0.8665 1.2908 0.0000 1.1197
#> 
#> Longitudinal outcome: y (family = beta, link = logit)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.3458 0.1107 -2.5599 -2.1296 0.0000 1.0026
#> sexfemale      -0.0676 0.1515 -0.3640  0.2273 0.6604 1.0060
#> time            0.3275 0.0521  0.2260  0.4297 0.0000 1.0017
#> sexfemale:time -0.0618 0.0738 -0.2060  0.0825 0.3953 1.0043
#> sigma           6.3279 0.4085  5.5334  7.1724 0.0000 1.0184
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 18 sec
```

[Back to top](#top)

### Censored linear mixed models

Some continuous longitudinal outcomes may have a censored nature. A
typical example of such outcomes is when we have a limit of detection
issue. That is, the values of the outcome cannot be detected below a
specified threshold having to do with the (laboratory) equipment used to
determine the measurements. In these settings, even if the complete data
follows a normal distribution the observed censored data cannot be
analyzed using a standard mixed model. The
[`mixed_model()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md)
function can accommodate such outcomes using the
[`censored.normal()`](https://rdrr.io/pkg/GLMMadaptive/man/extra_fams.html)
family object.

The following code simulates data from a joint model with a linear mixed
model for the longitudinal outcomes but applies censoring in the
realized longitudinal observations.

``` r

set.seed(1234)
n <- 200 # number of subjects
K <- 12 # number of measurements per subject
t_max <- 14 # maximum follow-up time

# we construct a data frame with the design:
# everyone has a baseline measurement, and then measurements at random 
# follow-up times up to t_max
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ time, data = DF)

betas <- c(-2.2, -0.25, 0.24, -0.05) # fixed effects coefficients
sigma <- 0.5 # errors' standard deviation
D11 <- 1.0 # variance of random intercepts
D22 <- 0.5 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
# we simulate normal longitudinal data
DF$y <- rnorm(n * K, mean = eta_y, sd = sigma)
# we assume that values below -4 are not observed, and set equal to -4
DF$ind <- as.numeric(DF$y < -4)
DF$y <- pmax(DF$y, -4)

upp_Cens <- 15 # fixed Type I censoring time
shape_wb <- 5 # shape Weibull
alpha <- 0.8 # association coefficients
gammas <- c("(Intercept)" = -9, "sex" = 0.5)
W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
# linear predictor for the survival model
eta_t <- as.vector(W %*% gammas)
# to simulate event times we use inverse transform sampling
# (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
# to find t, such that S(t) = u, where S(.) is the survival function, and u a 
# number from the Unif(0, 1) distribution. The function below calculates 
# log(u) - log(S(t)), and for a given u, we want to find t for which it equals
# zero. We do that below using the uniroot() function
invS <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h <- function (s) {
    X_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f <- as.vector(X_at_s %*% betas +
                     rowSums(Z_at_s * b[rep(i, nrow(Z_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f * alpha)
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h, lower = 0, upper = t)$value + log(u[i])
}
# we simulate the event times
u <- runif(n)
trueTimes <- numeric(n)
for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}

# we use fixed Type I right censoring denoting the end of the trial.
Ctimes <- upp_Cens
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator

# we keep the longitudinal measurements before the event times
DF$Time <- Time[DF$id]
DF$event <- event[DF$id]
DF <- DF[DF$time <= DF$Time, ]
```

The corresponding joint model is fitted with the following syntax:

``` r

DF_id <- DF[!duplicated(DF$id), ]
Cox_fit <- coxph(Surv(Time, event) ~ sex, data = DF_id)
CensNorm_MixMod <-
    mixed_model(cbind(y, ind) ~ sex * time, random = ~ time | id, data = DF,
                family = censored.normal())

jointFit <- jm(Cox_fit, CensNorm_MixMod, time_var = "time")
summary(jointFit)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = Cox_fit, Mixed_objects = CensNorm_MixMod, 
#>     time_var = "time")
#> 
#> Data Descriptives:
#> Number of groups: 200        Number of events: 165 (82.5%)
#> Number of observations:
#>   cbind(y, ind): 1346
#> 
#>                  DIC     WAIC      LPML
#> marginal    3721.930 4649.925 -3253.137
#> conditional 3335.145 3221.196 -1750.729
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.9504 (Intr)
#> time   0.6675 0.1453
#> 
#> Survival outcome:
#>                        Mean  StDev   2.5%  97.5%      P   Rhat
#> sexfemale            0.6167 0.2474 0.1442 1.1087 0.0111 1.0127
#> value(cbind(y, ind)) 0.8778 0.0659 0.7586 1.0214 0.0000 1.0444
#> 
#> Longitudinal outcome: cbind(y, ind) (family = censored normal, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.2588 0.1015 -2.4564 -2.0582 0.0000 1.0015
#> sexfemale      -0.0284 0.1440 -0.3145  0.2530 0.8489 1.0029
#> time            0.3335 0.0698  0.1973  0.4688 0.0000 1.0007
#> sexfemale:time -0.1174 0.0980 -0.3117  0.0692 0.2253 1.0012
#> sigma           0.4948 0.0137  0.4687  0.5215 0.0000 1.0020
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 14 sec
```

### Students’s-t mixed models

Outlying observations are a common issue in practice. Several methods
have been proposed in the literature for identifying such observations
in the context of longitudinal data. However, removing such values from
the analysis is generally not recommended unless we also have external
information as to why these values are outlying. Hence, we would need to
fit mixed models to accommodate these observations in these settings. A
well-known approach to achieve this is replacing the normal distribution
for the error terms in the linear mixed model with a Student’s-t
distribution with heavier tails.

The following syntax simulates data from a joint model with a
Student’s-t mixed effects model:

``` r

set.seed(1234)
n <- 200 # number of subjects
K <- 12 # number of measurements per subject
t_max <- 14 # maximum follow-up time

# we construct a data frame with the design:
# everyone has a baseline measurement, and then measurements at random 
# follow-up times up to t_max
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ time, data = DF)

betas <- c(-2.2, -0.25, 0.24, -0.05) # fixed effects coefficients
sigma <- 0.5 # error standard deviation
D11 <- 1.0 # variance of random intercepts
D22 <- 0.5 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
# we simulate Student's-t longitudinal data
DF$y <- eta_y + sigma * rt(n * K, df = 4)

upp_Cens <- 15 # fixed Type I censoring time
shape_wb <- 5 # shape Weibull
alpha <- 0.8 # association coefficients
gammas <- c("(Intercept)" = -9, "sex" = 0.5)
W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
# linear predictor for the survival model
eta_t <- as.vector(W %*% gammas)
# to simulate event times we use inverse transform sampling
# (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
# to find t, such that S(t) = u, where S(.) is the survival function, and u a 
# number from the Unif(0, 1) distribution. The function below calculates 
# log(u) - log(S(t)), and for a given u, we want to find t for which it equals
# zero. We do that below using the uniroot() function
invS <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h <- function (s) {
    X_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f <- as.vector(X_at_s %*% betas +
                     rowSums(Z_at_s * b[rep(i, nrow(Z_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f * alpha)
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h, lower = 0, upper = t)$value + log(u[i])
}
# we simulate the event times
u <- runif(n)
trueTimes <- numeric(n)
for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}

# we use fixed Type I right censoring denoting the end of the trial.
Ctimes <- upp_Cens
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator

# we keep the longitudinal measurements before the event times
DF$Time <- Time[DF$id]
DF$event <- event[DF$id]
DF <- DF[DF$time <= DF$Time, ]
```

To fit the corresponding joint model we use the
[`students.t()`](https://rdrr.io/pkg/GLMMadaptive/man/extra_fams.html)
family object in the call to
[`mixed_model()`](https://drizopoulos.github.io/JMbayes2/reference/sliced_model_generics.md):

``` r

DF_id <- DF[!duplicated(DF$id), ]
Cox_fit <- coxph(Surv(Time, event) ~ sex, data = DF_id)
Stdt_MixMod <-
    mixed_model(y ~ sex * time, random = ~ time | id, data = DF,
                family = students.t(df = 4))

jointFit <- jm(Cox_fit, Stdt_MixMod, time_var = "time")
summary(jointFit)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = Cox_fit, Mixed_objects = Stdt_MixMod, 
#>     time_var = "time")
#> 
#> Data Descriptives:
#> Number of groups: 200        Number of events: 165 (82.5%)
#> Number of observations:
#>   y: 1347
#> 
#>                  DIC     WAIC      LPML
#> marginal    5204.080 6855.990 -3958.089
#> conditional 4452.636 4352.831 -2347.894
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.9786 (Intr)
#> time   0.6777 0.1152
#> 
#> Survival outcome:
#>             Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale 0.0158 0.2483 -0.4729 0.4872 0.9431 1.0007
#> value(y)  0.8263 0.0704  0.6987 0.9883 0.0000 1.0382
#> 
#> Longitudinal outcome: y (family = Student's-t, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.2450 0.1073 -2.4526 -2.0341 0.0000 1.0117
#> sexfemale      -0.1222 0.1515 -0.4207  0.1754 0.4220 1.0056
#> time            0.2803 0.0711  0.1440  0.4204 0.0002 1.0035
#> sexfemale:time -0.0317 0.0995 -0.2268  0.1582 0.7460 1.0017
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 13 sec
```

[Back to top](#top)

### Negative binomial mixed models

Count longitudinal outcomes are typically modeled with the Poisson
distribution. However, these outcomes often exhibit more variance than
what is allowed from the Poisson distribution, leading to the well-known
problem of over-dispersion. To accommodate this over-dispersion,
typically, the negative binomial distribution is used.

The following piece of code simulates data from a joint model for count
longitudinal data that follow the negative binomial distribution:

``` r

set.seed(1234)
n <- 500 # number of subjects
K <- 10 # number of measurements per subject
t_max <- 5 # maximum follow-up time

# we construct a data frame with the design:
# everyone has a baseline measurement, and then measurements at random 
# follow-up times up to t_max
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ time, data = DF)

betas <- c(0.8, -0.5, 0.8, -0.5) # fixed effects coefficients
shape <- 2 # shape/size parameter of the negative binomial distribution
D11 <- 1.0 # variance of random intercepts
D22 <- 0.3 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
# mean of the Beta distribution
mu_y <- plogis(eta_y) # plogis(eta_y) = exp(eta_y) / (1 + exp(eta_y))
# we simulate negative binomial longitudinal data
DF$y <- rnbinom(n * K, size = shape, mu = exp(eta_y))

# simulate event times
upp_Cens <- 5 # fixed Type I censoring time
shape_wb <- 5 # shape Weibull
alpha <- 0.8 # association coefficient
gammas <- c("(Intercept)" = -9, "sex" = 0.5)
W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
# linear predictor for the survival model
eta_t <- as.vector(W %*% gammas)
# to simulate event times we use inverse transform sampling
# (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
# to find t, such that S(t) = u, where S(.) is the survival function, and u a 
# number from the Unif(0, 1) distribution. The function below calculates 
# log(u) - log(S(t)), and for a given u, we want to find t for which it equals
# zero. We do that below using the uniroot() function
invS <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h <- function (s) {
    X_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f <- as.vector(X_at_s %*% betas +
                     rowSums(Z_at_s * b[rep(i, nrow(Z_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f * alpha)
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h, lower = 0, upper = t)$value + log(u[i])
}
# we simulate the event times
u <- runif(n)
trueTimes <- numeric(n)
for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}

# we use fixed Type I right censoring denoting the end of the trial.
Ctimes <- upp_Cens
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator

# we keep the longitudinal measurements before the event times
DF$Time <- Time[DF$id]
DF$event <- event[DF$id]
DF <- DF[DF$time <= DF$Time, ]
```

The corresponding joint model is the fitted using the following syntax:

``` r

DF_id <- DF[!duplicated(DF$id), ]
Cox_fit <- coxph(Surv(Time, event) ~ sex, data = DF_id)
NB_MixMod <- mixed_model(y ~ sex * time, random = ~ time | id, data = DF,
                         family = GLMMadaptive::negative.binomial())

jointFit <- jm(Cox_fit, NB_MixMod, time_var = "time")
summary(jointFit)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = Cox_fit, Mixed_objects = NB_MixMod, 
#>     time_var = "time")
#> 
#> Data Descriptives:
#> Number of groups: 500        Number of events: 409 (81.8%)
#> Number of observations:
#>   y: 3842
#> 
#>                  DIC     WAIC      LPML
#> marginal    21350.36 21302.96 -10651.62
#> conditional 22646.19 22354.90 -11545.31
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 1.0498 (Intr) 
#> time   0.5402 -0.0625
#> 
#> Survival outcome:
#>             Mean  StDev   2.5%  97.5% P   Rhat
#> sexfemale 0.6153 0.1637 0.3055 0.9401 0 1.0198
#> value(y)  0.8136 0.0536 0.7180 0.9317 0 1.0447
#> 
#> Longitudinal outcome: y (family = negative binomial, link = log)
#>                   Mean  StDev    2.5%   97.5% P   Rhat
#> (Intercept)     0.8409 0.0799  0.6837  0.9953 0 1.0042
#> sexfemale      -0.5854 0.1144 -0.8088 -0.3657 0 1.0029
#> time            0.8621 0.0420  0.7812  0.9447 0 1.0057
#> sexfemale:time -0.5143 0.0584 -0.6286 -0.4016 0 1.0003
#> sigma           1.9916 0.0830  1.8320  2.1573 0 1.0020
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 29 sec
```

[Back to top](#top)

### Beta-binomial longitudinal outcomes

For count data and binomial data, we may have an over-dispersion
problem. To accommodate this, we can change the standard binomial
distribution to a beta-binomial one.

The following piece of code simulates data from a joint model for
binomial longitudinal data that follow the beta-binomial distribution:

``` r

set.seed(1234)
n <- 500 # number of subjects
K <- 8 # number of measurements per subject
t_max <- 10 # maximum follow-up time

# we construct a data frame with the design:
# everyone has a baseline measurement, and then measurements at random 
# follow-up times up to t_max
DF <- data.frame(id = rep(seq_len(n), each = K),
                 time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                 sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))

# design matrices for the fixed and random effects
X <- model.matrix(~ sex * time, data = DF)
Z <- model.matrix(~ time, data = DF)

betas <- c(-2.2, -0.25, 0.24, -0.05) # fixed effects coefficients
phi <- 5 # precision parameter of the Beta distribution
D11 <- 1.0 # variance of random intercepts
D22 <- 0.5 # variance of random slopes

# we simulate random effects
b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)))
# linear predictor
eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, ]))
# mean of the Beta distribution
mu_y <- plogis(eta_y) # plogis(eta_y) = exp(eta_y) / (1 + exp(eta_y))
# we simulate probabilities from the Beta distribution
probs <- rbeta(n * K, shape1 = mu_y * phi, shape2 = phi * (1 - mu_y))
# we transform to (0, 1)
probs <- (probs * (nrow(DF) - 1) + 0.5) / nrow(DF)
# we simulate binomial data use the probs
DF$y <- rbinom(n * K, size = 20, prob = probs)

upp_Cens <- 15 # fixed Type I censoring time
shape_wb <- 5 # shape Weibull
alpha <- 0.8 # association coefficients
gammas <- c("(Intercept)" = -9, "sex" = 0.5)
W <- model.matrix(~ sex, data = DF[!duplicated(DF$id), ])
# linear predictor for the survival model
eta_t <- as.vector(W %*% gammas)
# to simulate event times we use inverse transform sampling
# (https://en.wikipedia.org/wiki/Inverse_transform_sampling). Namely, we want 
# to find t, such that S(t) = u, where S(.) is the survival function, and u a 
# number from the Unif(0, 1) distribution. The function below calculates 
# log(u) - log(S(t)), and for a given u, we want to find t for which it equals
# zero. We do that below using the uniroot() function
invS <- function (t, i) {
  # i denotes the subject
  sex_i <- W[i, 2L]
  # h() is the hazard function and we assume a Weibull baseline hazard
  h <- function (s) {
    X_at_s <- cbind(1, sex_i, s, sex_i * s)
    Z_at_s <- cbind(1, s)
    # the linear predictor from the mixed model evaluated at time s
    f <- as.vector(X_at_s %*% betas +
                     rowSums(Z_at_s * b[rep(i, nrow(Z_at_s)), ]))
    exp(log(shape_wb) + (shape_wb - 1) * log(s) + eta_t[i] + f * alpha)
  }
  # -log(S(t)) = H(t), where H(t) is the cumulative hazard function
  integrate(h, lower = 0, upper = t)$value + log(u[i])
}
# we simulate the event times
u <- runif(n)
trueTimes <- numeric(n)
for (i in seq_len(n)) {
    Up <- 100
    Root <- try(uniroot(invS, interval = c(1e-05, Up), i = i)$root, TRUE)
    trueTimes[i] <- if (!inherits(Root, "try-error")) Root else 150
}

# we use fixed Type I right censoring denoting the end of the trial.
Ctimes <- upp_Cens
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator

# we keep the longitudinal measurements before the event times
DF$Time <- Time[DF$id]
DF$event <- event[DF$id]
DF <- DF[DF$time <= DF$Time, ]
```

The corresponding joint model is then fitted with the syntax:

``` r

DF_id <- DF[!duplicated(DF$id), ]
Cox_fit <- coxph(Surv(Time, event) ~ sex, data = DF_id)
BetaBinom_MixMod <-
    mixed_model(cbind(y, 20 - y) ~ sex * time, random = ~ time | id, data = DF,
                family = beta.binomial())

jointFit <- jm(Cox_fit, BetaBinom_MixMod, time_var = "time")
summary(jointFit)
#> 
#> Call:
#> JMbayes2::jm(Surv_object = Cox_fit, Mixed_objects = BetaBinom_MixMod, 
#>     time_var = "time")
#> 
#> Data Descriptives:
#> Number of groups: 500        Number of events: 395 (79%)
#> Number of observations:
#>   cbind(y, 20 - y): 2837
#> 
#>                  DIC     WAIC      LPML
#> marginal    11924.23 11874.32 -5961.406
#> conditional 13497.76 13280.00 -7050.638
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 1.0211 (Intr) 
#> time   0.7168 -0.0097
#> 
#> Survival outcome:
#>                           Mean StDev   2.5%  97.5%     P   Rhat
#> sexfemale               0.5831 0.194 0.2189 0.9635 2e-04 1.0100
#> value(cbind(y, 20 - y)) 0.9353 0.069 0.8212 1.0879 0e+00 1.1183
#> 
#> Longitudinal outcome: cbind(y, 20 - y) (family = beta binomial, link = logit)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.1601 0.0958 -2.3489 -1.9765 0.0000 1.0080
#> sexfemale      -0.2744 0.1311 -0.5363 -0.0163 0.0378 1.0075
#> time            0.3313 0.0524  0.2291  0.4335 0.0000 1.0025
#> sexfemale:time -0.1793 0.0736 -0.3258 -0.0353 0.0140 1.0016
#> sigma           4.8418 0.2774  4.3026  5.3913 0.0000 1.0440
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 30 sec
```

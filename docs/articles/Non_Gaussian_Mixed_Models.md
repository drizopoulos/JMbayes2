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
[`beta.fam()`](https://drizopoulos.github.io/GLMMadaptive/reference/extra_fams.html)
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
#> marginal    -3700.436 -3722.514 1751.861
#> conditional -4171.526 -3943.855 1819.955
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.8443 (Intr)
#> time   0.4663 0.0785
#> 
#> Survival outcome:
#>             Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale 0.2169 0.2916 -0.3701 0.7895 0.4604 1.0055
#> value(y)  1.0584 0.0948  0.8890 1.2615 0.0000 1.0643
#> 
#> Longitudinal outcome: y (family = beta, link = logit)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.3446 0.1085 -2.5563 -2.1344 0.0000 1.0111
#> sexfemale      -0.0625 0.1511 -0.3630  0.2353 0.6787 1.0099
#> time            0.3272 0.0515  0.2271  0.4287 0.0000 1.0062
#> sexfemale:time -0.0631 0.0733 -0.2072  0.0804 0.3842 1.0084
#> sigma           6.2935 0.3824  5.5669  7.0915 0.0000 1.0315
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 17 sec
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
[`censored.normal()`](https://drizopoulos.github.io/GLMMadaptive/reference/extra_fams.html)
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
#> marginal    4958.446 7521.087 -4660.703
#> conditional 2273.742 3219.776 -1756.184
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.9478 (Intr)
#> time   0.6600 0.1456
#> 
#> Survival outcome:
#>                        Mean  StDev   2.5% 97.5%      P   Rhat
#> sexfemale            0.5991 0.2473 0.1252 1.073 0.0142 1.0043
#> value(cbind(y, ind)) 0.8672 0.0657 0.7457 1.007 0.0000 1.0359
#> 
#> Longitudinal outcome: cbind(y, ind) (family = censored normal, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.2608 0.1012 -2.4567 -2.0624 0.0000 1.0016
#> sexfemale      -0.0255 0.1438 -0.3087  0.2568 0.8616 1.0020
#> time            0.3334 0.0690  0.1988  0.4679 0.0000 1.0004
#> sexfemale:time -0.1179 0.0973 -0.3099  0.0708 0.2307 1.0003
#> sigma           0.4949 0.0138  0.4692  0.5230 0.0000 1.0003
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
[`students.t()`](https://drizopoulos.github.io/GLMMadaptive/reference/extra_fams.html)
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
#> marginal    6151.498 8249.053 -4370.657
#> conditional 3587.961 4350.715 -2338.512
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.9795 (Intr)
#> time   0.6721 0.1166
#> 
#> Survival outcome:
#>            Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale 0.017 0.2466 -0.4614 0.5010 0.9307 1.0058
#> value(y)  0.823 0.0645  0.7074 0.9625 0.0000 1.0910
#> 
#> Longitudinal outcome: y (family = Student's-t, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.2460 0.1074 -2.4519 -2.0340 0.0000 1.0057
#> sexfemale      -0.1217 0.1511 -0.4158  0.1718 0.4153 1.0029
#> time            0.2794 0.0705  0.1430  0.4168 0.0002 1.0066
#> sexfemale:time -0.0312 0.0985 -0.2256  0.1595 0.7489 1.0049
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
#> marginal    21716.78 21678.49 -10839.52
#> conditional 22255.35 22347.60 -11559.46
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 1.0398 (Intr) 
#> time   0.5364 -0.0502
#> 
#> Survival outcome:
#>             Mean  StDev   2.5%  97.5% P   Rhat
#> sexfemale 0.6168 0.1653 0.3061 0.9410 0 1.0150
#> value(y)  0.8175 0.0542 0.7129 0.9267 0 1.0214
#> 
#> Longitudinal outcome: y (family = negative binomial, link = log)
#>                   Mean  StDev    2.5%   97.5% P   Rhat
#> (Intercept)     0.8415 0.0785  0.6851  0.9938 0 1.0033
#> sexfemale      -0.5827 0.1134 -0.8076 -0.3660 0 1.0009
#> time            0.8620 0.0414  0.7813  0.9428 0 1.0017
#> sexfemale:time -0.5149 0.0582 -0.6302 -0.4018 0 1.0008
#> sigma           1.9878 0.0837  1.8274  2.1574 0 1.0021
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
#>                   DIC     WAIC       LPML
#> marginal    16068.827 21506.46 -10946.349
#> conditional  9850.636 13256.09  -7023.212
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.9996 (Intr)
#> time   0.7069 0.0040
#> 
#> Survival outcome:
#>                           Mean  StDev   2.5%  97.5%      P   Rhat
#> sexfemale               0.5880 0.1985 0.2080 0.9867 0.0053 1.0093
#> value(cbind(y, 20 - y)) 0.9487 0.0722 0.8227 1.1051 0.0000 1.2894
#> 
#> Longitudinal outcome: cbind(y, 20 - y) (family = beta binomial, link = logit)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.1503 0.0946 -2.3383 -1.9667 0.0000 1.0373
#> sexfemale      -0.2702 0.1313 -0.5340 -0.0161 0.0373 1.0035
#> time            0.3288 0.0524  0.2276  0.4298 0.0000 1.0048
#> sexfemale:time -0.1778 0.0742 -0.3239 -0.0336 0.0187 1.0009
#> sigma           4.7852 0.2768  4.2656  5.3463 0.0000 1.0598
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 30 sec
```

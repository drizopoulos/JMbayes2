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
#> Number of Groups: 200        Number of events: 158 (79%)
#> Number of Observations:
#>   y: 1182
#> 
#>                   DIC      WAIC     LPML
#> marginal    -4018.171 -4061.024 1910.317
#> conditional -3837.701 -3938.538 1799.508
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.8436 (Intr)
#> time   0.4664 0.0812
#> 
#> Survival Outcome:
#>             Mean  StDev    2.5%  97.5%      P   Rhat
#> sexfemale 0.2153 0.2944 -0.3672 0.7865 0.4669 1.0030
#> value(y)  1.0656 0.1018  0.8866 1.2781 0.0000 1.0832
#> 
#> Longitudinal Outcome: y (family = beta, link = logit)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.3410 0.1115 -2.5592 -2.1270 0.0000 1.0046
#> sexfemale      -0.0689 0.1486 -0.3584  0.2241 0.6462 1.0144
#> time            0.3271 0.0523  0.2258  0.4312 0.0000 1.0022
#> sexfemale:time -0.0615 0.0729 -0.2035  0.0801 0.3898 1.0059
#> sigma           6.2757 0.3932  5.5169  7.0813 0.0000 1.0064
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 20 sec
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
#> Number of Groups: 200        Number of events: 165 (82.5%)
#> Number of Observations:
#>   cbind(y, ind): 1346
#> 
#>                  DIC     WAIC      LPML
#> marginal    3718.336 4647.791 -3267.595
#> conditional 3330.882 3218.840 -1744.693
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.9502 (Intr)
#> time   0.6633 0.1447
#> 
#> Survival Outcome:
#>                        Mean  StDev   2.5%  97.5%      P   Rhat
#> sexfemale            0.6065 0.2483 0.1192 1.0889 0.0107 1.0027
#> value(cbind(y, ind)) 0.8702 0.0641 0.7540 1.0005 0.0000 1.0477
#> 
#> Longitudinal Outcome: cbind(y, ind) (family = censored normal, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.2588 0.1023 -2.4597 -2.0543 0.0000 1.0019
#> sexfemale      -0.0300 0.1433 -0.3120  0.2483 0.8382 1.0011
#> time            0.3344 0.0698  0.1989  0.4710 0.0000 1.0001
#> sexfemale:time -0.1175 0.0983 -0.3082  0.0732 0.2260 1.0011
#> sigma           0.4947 0.0138  0.4690  0.5227 0.0000 1.0028
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 15 sec
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
#> Number of Groups: 200        Number of events: 165 (82.5%)
#> Number of Observations:
#>   y: 1347
#> 
#>                  DIC     WAIC      LPML
#> marginal    5197.868 6842.793 -3988.557
#> conditional 4447.164 4349.718 -2333.522
#> 
#> Random-effects covariance matrix:
#>                     
#>        StdDev   Corr
#> (Intr) 0.9791 (Intr)
#> time   0.6749 0.1175
#> 
#> Survival Outcome:
#>             Mean  StDev    2.5%  97.5%     P   Rhat
#> sexfemale 0.0207 0.2538 -0.4730 0.5131 0.932 1.0012
#> value(y)  0.8350 0.0647  0.7152 0.9658 0.000 1.0569
#> 
#> Longitudinal Outcome: y (family = Student's-t, link = identity)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.2485 0.1083 -2.4626 -2.0334 0.0000 1.0038
#> sexfemale      -0.1212 0.1513 -0.4205  0.1691 0.4267 1.0005
#> time            0.2828 0.0710  0.1455  0.4235 0.0000 1.0007
#> sexfemale:time -0.0335 0.0993 -0.2263  0.1598 0.7318 1.0007
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 15 sec
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
#> Number of Groups: 500        Number of events: 409 (81.8%)
#> Number of Observations:
#>   y: 3842
#> 
#>                  DIC     WAIC      LPML
#> marginal    21347.40 21300.55 -10650.44
#> conditional 22639.19 22350.27 -11565.82
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 1.0401 (Intr) 
#> time   0.5382 -0.0551
#> 
#> Survival Outcome:
#>             Mean  StDev   2.5%  97.5% P   Rhat
#> sexfemale 0.6141 0.1608 0.2967 0.9385 0 1.0126
#> value(y)  0.8114 0.0536 0.7136 0.9274 0 1.0458
#> 
#> Longitudinal Outcome: y (family = negative binomial, link = log)
#>                   Mean  StDev    2.5%   97.5% P   Rhat
#> (Intercept)     0.8422 0.0791  0.6852  0.9952 0 1.0116
#> sexfemale      -0.5812 0.1143 -0.8014 -0.3573 0 1.0035
#> time            0.8617 0.0415  0.7813  0.9427 0 1.0053
#> sexfemale:time -0.5161 0.0591 -0.6298 -0.3989 0 1.0010
#> sigma           1.9921 0.0848  1.8281  2.1608 0 1.0027
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 27 sec
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
#> Number of Groups: 500        Number of events: 395 (79%)
#> Number of Observations:
#>   cbind(y, 20 - y): 2837
#> 
#>                  DIC     WAIC      LPML
#> marginal    11921.66 11871.70 -5942.471
#> conditional 13491.02 13297.62 -7547.969
#> 
#> Random-effects covariance matrix:
#>                      
#>        StdDev   Corr 
#> (Intr) 1.0077 (Intr) 
#> time   0.7146 -0.0045
#> 
#> Survival Outcome:
#>                           Mean  StDev   2.5%  97.5%      P   Rhat
#> sexfemale               0.5800 0.1952 0.2149 0.9642 0.0036 1.0050
#> value(cbind(y, 20 - y)) 0.9346 0.0598 0.8248 1.0519 0.0000 1.0299
#> 
#> Longitudinal Outcome: cbind(y, 20 - y) (family = beta binomial, link = logit)
#>                   Mean  StDev    2.5%   97.5%      P   Rhat
#> (Intercept)    -2.1593 0.0940 -2.3486 -1.9754 0.0000 1.0194
#> sexfemale      -0.2722 0.1310 -0.5360 -0.0178 0.0336 1.0087
#> time            0.3316 0.0522  0.2299  0.4334 0.0000 1.0039
#> sexfemale:time -0.1798 0.0749 -0.3276 -0.0368 0.0142 1.0007
#> sigma           4.8145 0.2810  4.2778  5.3776 0.0000 1.0262
#> 
#> MCMC summary:
#> chains: 3 
#> iterations per chain: 3500 
#> burn-in per chain: 500 
#> thinning: 1 
#> time: 50 sec
```

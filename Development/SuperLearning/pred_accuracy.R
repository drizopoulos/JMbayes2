library("JMbayes2")
library("mice")

# Function to simulate from a joint model
simulate_jm <- function (n, seed) {
    # random seed
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1L)
    RNGstate <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
    set.seed(seed)
    K <- 12 # number of measurements per subject
    t_max <- 14 # maximum follow-up time

    # we construct a data frame with the design:
    # everyone has a baseline measurement, and then measurements at random
    # follow-up times up to t_max
    DF <- data.frame(id = rep(seq_len(n), each = K),
                     # time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
                     time = rep(seq(0, t_max, length.out = K), n),
                     treat = rep(gl(2, n/2, labels = c("placebo", "active")), each = K))

    # design matrices for the fixed and random effects
    X <- model.matrix(~ treat * time, data = DF)
    Z <- model.matrix(~ time, data = DF)

    betas <- c(-2.2, -0.25, 0.24, -0.45) # fixed effects coefficients
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

    upp_Cens <- 15 # fixed Type I censoring time
    shape_wb <- 5 # shape Weibull
    alpha <- 0.8 # association coefficients
    gammas <- c("(Intercept)" = -9, "treat" = 0.8)
    W <- model.matrix(~ treat, data = DF[!duplicated(DF$id), ])
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
    DF[DF$time <= DF$Time, ]
}

# Simulate some training and some test data
training <- simulate_jm(n = 200L, seed = 1234L)
test <- simulate_jm(n = 200L, seed = 4321L)

# Fit the joint model
training_id <- training[!duplicated(training$id), ]
Cox_fit <- coxph(Surv(Time, event) ~ treat, data = training_id)
lme_fit <- lme(y ~ treat * time, data = training, random = ~ time | id)
joint_fit <- jm(Cox_fit, lme_fit, time_var = "time")
summary(joint_fit)

# We want to calculate causal predictions for subjects at risk at time t0 in
# the test data. We only keep the measurements before t0
t0 <- 5
aTest <- test[test$Time > t0 & test$time <= t0, ]

# We turn into a wide dataset
aTest_wide <- reshape(aTest, timevar = "time", direction = "wide",
                      idvar = "id", v.names = "y", sep = "_")

# We turn event into a factor
aTest_wide$event <- factor(aTest_wide$event, levels = 0:1,
                           labels = c('alive', 'dead'))

# create a copy of aTest_wide with switch treatment and no event time
# and event status
copy_aTest_wide <- aTest_wide
lvs_treat <- levels(aTest_wide$treat)
copy_aTest_wide$treat <- ifelse(copy_aTest_wide$treat == lvs_treat[1L],
                                lvs_treat[2L], lvs_treat[1L])
# We want to impute the event times and event status under the counterfactual
# scenario of opposite treatments. Hence, we set these variables to NA
copy_aTest_wide$event <- NA
copy_aTest_wide$Time <- NA

#############################################################################

# impute event times and event information
imps <- mice(data = rbind(aTest_wide[-1], copy_aTest_wide[-1]), maxit = 0L)
impsmeth <- imps$meth
impsmeth[c('Time', 'event')] <- c('rf', 'logreg')
imps <- mice(data = rbind(aTest_wide[-1], copy_aTest_wide[-1]), seed = 2023L,
             printFlag = FALSE, methods = impsmeth)


D_i <- complete(imps)
D_i$syn <- rep(0:1, each = nrow(aTest_wide))
D_i$id <- c(aTest_wide$id, paste0(aTest_wide$id, "_clone"))

vnams <- names(D_i)[grep("y_", names(D_i), fixed = TRUE)]
longD_i <- reshape(D_i, direction = "long", varying = vnams, sep = "_")



################################################################################
################################################################################
################################################################################


# We will work with discrete times. We create a time grid of length 0.1
times_grid <- seq(0, max(aTest$Time) + 1, by = 0.1)
# We extract per subject the times longitudinal measurements were taken
times <- with(aTest, split(time, id))
Time <- with(aTest, split(Time, id))

# We find in which interval in the time grid each time belongs
mapply2 <- function (...) mapply(..., SIMPLIFY = FALSE)
aTest$time_int <-
    unlist(mapply2(findInterval, x = times, MoreArgs = list(vec = times_grid)))
aTest$Time_int <-
    unlist(mapply(findInterval, x = Time, MoreArgs = list(vec = times_grid)))
# We turn 'Time_int' into a ordered factor, such that to use the
# proportional odds model for the imputation
aTest$Time_int <- ordered(aTest$Time_int)
# We also turn event into a factor
aTest$event <- factor(aTest$event, levels = 0:1, labels = c('alive', 'dead'))
# we keep the variables we will need
aTest <- aTest[c("id", "time_int", "y", "Time_int", "event", "treat")]
# We turn into a wide dataset (this gives a warning for multiple matches for
# specific time points, i.e., for some subjects we have more than one
# longitudinal measurements in a time interval specified in 'times_grid'. Now
# only the first value is taken. This is to be fixed, e.g., calculating the
# average of the longitudinal measurements in the interval)
aTest_wide <- reshape(aTest, timevar = "time_int", direction = "wide",
                      idvar = "id", v.names = "y", sep = "_")
# We order the columns per time (just for cosmetic reasons)
all_cols <- paste0("y_", 1:200)
aTest_wide <- aTest_wide[c("id", "Time_int", "event", "treat",
                           all_cols[all_cols %in% names(aTest_wide)])]

# create a copy of aTest_wide with switch treatment and no event time
# and event status
copy_aTest_wide <- aTest_wide
lvs_treat <- levels(aTest_wide_wide$treat)
copy_aTest_wide$treat[copy_aTest_wide$treat == lvs_treat[1L]] <- lvs_treat[2L]
copy_aTest_wide$treat[copy_aTest_wide$treat == lvs_treat[2L]] <- lvs_treat[1L]
# We want to impute the event times and event status under the counterfactual
# scenario of opposite treatments. Hence, we set these variables to NA
copy_aTest_wide$event <- NA
copy_aTest_wide$Time_int <- NA


DF <- pbc2[c("id", "year", "serBilir", "years", "status", "drug", "sex", "age")]
times <- with(DF, split(year, id))
Time <- with(DF, split(years, id))
times_grid <- seq(0, 15, by = 0.5)
DF$times <- unlist(mapply(findInterval, x = times, MoreArgs = list(vec = times_grid)))
DF$Time <- unlist(mapply(findInterval, x = Time, MoreArgs = list(vec = times_grid)))
DF <- DF[c("id", "times", "serBilir", "Time", "status", "drug", "sex", "age")]

DF_wide <- reshape(DF, timevar = "times", direction = "wide", idvar = "id",
                   v.names = "serBilir", sep = "_")
all_cols <- paste0("serBilir_", 1:57)
DF_wide <- DF_wide[c("id", "Time", "status", "drug", "sex", "age",
                     all_cols[all_cols %in% names(DF_wide)])]

#############################################################################

# impute longitudinal data
imps_long <- mice(DF_wide[-1], seed = 2023L, printFlag = FALSE)
DD <- complete(imps_long)
View(DD)

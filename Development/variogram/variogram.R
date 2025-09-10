library("rbenchmark")

variogram1 <- function (y, times, id) {
    unq_id <- unique(id)
    n <- length(unq_id)
    vv <- vt <- vector("list", n)
    for (i in seq_len(n)) {
        id_i <- id == unq_id[i]
        y_i <- y[id_i]
        times_i <- times[id_i]
        n_i <- length(y_i)
        if (n_i > 1L) {
            inds <- combn(n_i, 2)
            vv[[i]] <- 0.5 * (y_i[inds[1L, ]] - y_i[inds[2L, ]])^2
            vt[[i]] <- abs(times_i[inds[1L, ]] - times_i[inds[2L, ]])
        } else {
            vv[[i]] <- vt[[i]] <- 0
        }
    }
    cbind(unlist(vt, use.names = FALSE), unlist(vv, use.names = FALSE))
}

variogram2 <- function (y, times, id) {
    variogram_cpp(split(y, id), split(times, id))
}

################################################################################
################################################################################

n <- 500
k <- 20
yy <- rnorm(n * k)
tt <- runif(n * k, 0, 150)
id <- rep(seq_len(n), each = k)

v1 <- variogram1(yy, tt, id)
v2 <- variogram2(yy, tt, id)

all.equal(v1, v2)


benchmark(R = variogram1(yy, tt, id),
          Cpp = variogram2(yy, tt, id),
          replications = 30)

data("pbc2", package = "JM")

v0 <- variogram(log(pbc2$serBilir), pbc2$year, unclass(pbc2$id), FALSE)[[1]]
v1 <- variogram1(log(pbc2$serBilir), pbc2$year, unclass(pbc2$id))
v2 <- variogram2(log(pbc2$serBilir), pbc2$year, unclass(pbc2$id))

all.equal(v0, v1, check.attributes = FALSE)
all.equal(v0, v2, check.attributes = FALSE)
all.equal(v1, v2)


benchmark(R = variogram1(log(pbc2$serBilir), pbc2$year, unclass(pbc2$id)),
          Cpp = variogram2(log(pbc2$serBilir), pbc2$year, unclass(pbc2$id)),
          replications = 500)



object = jointFit#jFit1_value
nsim = 40L
newdata = NULL
seed = 123L
process = "longitudinal"
type = "variogram"
outcomes = Inf
percentiles = c(0.025, 0.975)
random_effects = "mcmc"
Fforms_fun = NULL


out <- simulate(object, nsim = nsim, newdata = newdata,
             process = "longitudinal", include_outcome = TRUE,
             random_effects = random_effects, seed = seed,
             Fforms_fun = Fforms_fun)
yy <- out$outcome
out <- out[names(out) != "outcome"]
n_outcomes <- length(out)
index <- seq_len(n_outcomes)

j = 1
y <- yy[[j]]
X <- attr(y, "X")
lm_fit <- attr(y, "lm_fit")
resd_obs <- attr(y, "resd_obs")
tt <- attr(y, "times")
id <- attr(y, "id")


vrgm_DF <- variogram(resd_obs, tt, id)[[1]]
loess_obs <- loess(diffs2 ~ time_lag, data = data.frame(vrgm_DF))
ttt <- seq(min(vrgm_DF[, 1]), max(vrgm_DF[, 1]), len = 101)
vrgm_obs_loess <- predict(loess_obs, data.frame(time_lag = ttt))
vrgm_rep_loess <- matrix(0, length(ttt), ncol(out[[j]]))
for (i in seq_len(ncol(out[[j]]))) {
    na_ind <- is.na(out[[j]][, i])
    vrgm_DF <- variogram(out[[j]][!na_ind, i] -
                             lm_fit[!na_ind], tt[!na_ind], id[!na_ind])[[1L]]
    loess_rep_i <- loess(diffs2 ~ time_lag, data = data.frame(vrgm_DF))
    vrgm_rep_loess[, i] <- predict(loess_rep_i,
                                   data.frame(time_lag = ttt))
}
matplot(ttt, vrgm_rep_loess, type = "l", col = "lightgrey",
        lty = 1, xlab = "Time lags", ylab = "Half Squared Differences")
lines(ttt, vrgm_obs_loess, lwd = 2)


vrgm_DF <- variogram(resd_obs, tt, id)[[1]]
vrgm_obs_loess <-
    loess.smooth(vrgm_DF[, "time_lag"], vrgm_DF[, "diffs2"], degree = 2,
                 span = 0.75, family = "gaussian")
vrgm_rep_loess <- matrix(0, length(vrgm_obs_loess$y), ncol(out[[j]]))
for (i in seq_len(ncol(out[[j]]))) {
    not_na <- !is.na(out[[j]][, i])
    vrgm_DF <- variogram(out[[j]][not_na, i] -
                             lm_fit[not_na], tt[not_na], id[not_na])[[1L]]
    loess_rep_i <-
        loess.smooth(vrgm_DF[, "time_lag"], vrgm_DF[, "diffs2"], degree = 2,
                     span = 0.75, family = "gaussian")
    vrgm_rep_loess[, i] <- loess_rep_i$y
}
matplot(vrgm_obs_loess$x, vrgm_rep_loess, type = "l",
        col = "lightgrey", lty = 1,
        xlab = "Time lags", ylab = "Half Squared Differences")
lines(vrgm_obs_loess, lwd = 2)


xxx <- rnorm(500)
yyy <- (2 + 0.5 * xxx - 0.05 * xxx^2 + 0.08 * xxx^3 + rnorm(500))^2
loess_fit <- loess(yyy ~ xxx, data = data.frame(xxx = xxx, yyy = yyy),
                   span = 0.75)
fit1 <- predict(loess_fit, data.frame(xxx = seq(-3, 3, len = 201)))
fit2 <- loess.smooth(xxx, yyy, span = 0.75)$y
cbind(fit1, fit2)
plot(xxx, yyy)
ind <- order(xxx)
#lines(xxx[ind], fit1[ind], lwd = 2, col = "red")
lines(seq(-3, 3, len = 201), fit1, lwd = 2, col = "red")
lines(loess.smooth(xxx, yyy, span = 0.75, degree = 2), lwd = 2, col = "blue")




vrgm_DF <- variogram(resd_obs, tt, id)[[1]]
loess_obs <- loess(diffs2 ~ time_lag, data = data.frame(vrgm_DF), span = 0.75)
ttt <- seq(min(vrgm_DF[, 1]), max(vrgm_DF[, 1]), len = 101)
fit1 <- predict(loess_obs, data.frame(time_lag = ttt))

plot(diffs2 ~ time_lag, data = data.frame(vrgm_DF), pch = ".", ylim = c(0, 1))
lines(ttt, fit1, lwd = 2, col = "red")
lines(loess.smooth(vrgm_DF[, 1], vrgm_DF[, 2], span = 0.75, degree = 2, family = "g"),
      lwd = 2, col = "blue")





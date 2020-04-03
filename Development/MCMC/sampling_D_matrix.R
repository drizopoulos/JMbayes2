library("survival")
library("nlme")
library("GLMMadaptive")
library("splines")
data("pbc2", package = "JM")
data("pbc2.id", package = "JM")
source(file.path(getwd(), "R/jm.R"))
source(file.path(getwd(), "R/help_functions.R"))
source(file.path(getwd(), "Development/jm/R_to_Cpp.R"))
source(file.path(getwd(), "Development/jm/PBC_data.R"))

fm1 <- lme(log(serBilir) ~ year * sex + I(year^2) + age + prothrombin,
           data = pbc2, random = ~ year | id)
fm2 <- lme(serChol ~ ns(year, 3) + sex + age, data = pbc2, random = ~ year | id,
           na.action = na.exclude)
fm3 <- mixed_model(hepatomegaly ~ sex + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())
fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())
Mixed_objects <- list(fm1, fm2, fm3, fm4)

D_lis <- lapply(Mixed_objects, extract_D)
D <- bdiag(D_lis)
invD <- solve(D)

##########################################################################################
##########################################################################################
dmvnorm <- function (x, mu, Sigma = NULL, invSigma = NULL, log = TRUE,
          prop = TRUE) {
    if (!is.matrix(x))
        x <- rbind(x)
    if (is.matrix(mu)) {
        if (nrow(mu) != nrow(x))
            stop("incorrect dimensions for 'mu'.")
        p <- ncol(mu)
    }
    else {
        p <- length(mu)
        mu <- rep(mu, each = nrow(x))
    }
    if (is.null(Sigma) && is.null(invSigma))
        stop("'Sigma' or 'invSigma' must be given.")
    if (!is.null(Sigma)) {
        if (is.list(Sigma)) {
            ev <- Sigma$values
            evec <- Sigma$vectors
        } else {
            ed <- eigen(Sigma, symmetric = TRUE)
            ev <- ed$values
            evec <- ed$vectors
        }
        invSigma <- evec %*% (t(evec)/ev)
        if (!prop)
            logdetSigma <- sum(log(ev))
    } else {
        if (!prop)
            logdetSigma <- -determinant(as.matrix(invSigma))$modulus
    }
    ss <- x - mu
    quad <- 0.5 * rowSums((ss %*% invSigma) * ss)
    if (!prop)
        fact <- -0.5 * (p * log(2 * pi) + logdetSigma)
    if (log) {
        if (!prop)
            as.vector(fact - quad)
        else as.vector(-quad)
    } else {
        if (!prop)
            as.vector(exp(fact - quad))
        else as.vector(exp(-quad))
    }
}

cor2cov <- function (R, vars, sds = NULL) {
    p <- nrow(R)
    if (is.null(sds)) sds <- sqrt(vars)
    sds * R * rep(sds, each = p)
}

ddirichlet <- function (x, alpha, log = FALSE) {
    dirichlet1 <- function(x, alpha) {
        logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
        s <- sum((alpha - 1) * log(x))
        exp(sum(s) - logD)
    }
    if (!is.matrix(x))
        if (is.data.frame(x))
            x <- as.matrix(x)
        else x <- t(x)
        if (!is.matrix(alpha))
            alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x),
                            byrow = TRUE)
        if (any(dim(x) != dim(alpha)))
            stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
        pd <- vector(length = nrow(x))
        for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i,
                                                               ])
        pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
        pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
        if (log) return(log(pd)) else return(pd)
}

rdirichlet <- function (n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
    x / rowSums(x)
}

update_scale <- function (scale, rate, target_acc, it, c1 = 0.8, c0 = 1) {
    g1 <- (it + 1)^(1 - c1)
    g2 <- c0 * g1
    exp(log(scale) + g2 * (rate - target_acc))
}

dht <- function (x, sigma = 10, df = 1, log = FALSE) {
    ind <- x > 0
    out <- rep(as.numeric(NA), length(x))
    log_const <- log(2) + lgamma(0.5 * (df + 1)) - lgamma(0.5 * df) -
        0.5 * (log(df) + log(pi)) - log(sigma)
    log_kernel <- - 0.5 * (df + 1) * log(1 + x[ind]^2 / (df * sigma^2))
    out[ind] <- log_const + log_kernel
    if (log) out else exp(out)
}

##########################################################################################
##########################################################################################

# sampling using the global variance, simplex weights approach

p <- ncol(D)
R <- cov2cor(D)
vars <- diag(D)
tau <- sum(vars) / p # the trace is p * tau
simplex <- vars / sum(vars)

D <- cor2cov(R, p * tau * simplex)

b <- MASS::mvrnorm(1000, rep(0, p), D)

target_log_dist <- function (tau, simplex) {
    p <- length(simplex)
    D <- cor2cov(R, p * tau * simplex)
    log_p_b <- sum(dmvnorm(b, rep(0, p), D, log = TRUE, prop = FALSE))
    log_p_tau <- dgamma(tau, 1, 1, log = TRUE)
    log_p_simplex <- ddirichlet(simplex, rep(1, p), log = TRUE)
    log_p_b + log_p_tau + log_p_simplex
}

M <- 3000
acceptance_tau <- taus <- numeric(M)
acceptance_simplex <- numeric(M)
simplexes <- matrix(0.0, M, p)
current_tau <- tau
current_simplex <- simplex
scale_tau <- 0.04
scale_simplex <- 2e4
#scale_tau <- 0.1
#scale_simplex <- 1e3
for (m in seq_len(M)) {
    log_mu <- log(current_tau) - 0.5 * scale_tau^2
    proposed_tau <- rlnorm(1, log_mu, scale_tau)
    numerator <- target_log_dist(proposed_tau, current_simplex) +
        dlnorm(current_tau, log_mu, scale_tau, log = TRUE)
    denominator <- target_log_dist(current_tau, current_simplex) +
        dlnorm(proposed_tau, log_mu, scale_tau, log = TRUE)
    log_ratio <- numerator - denominator
    if (log_ratio > log(runif(1))) {
        current_tau <- proposed_tau
        acceptance_tau[m] <- 1
    }
    taus[m] <- current_tau
    #######
    proposed_simplex <- c(rdirichlet(1, scale_simplex * current_simplex))

    numerator <- target_log_dist(current_tau, proposed_simplex) +
        ddirichlet(current_simplex, current_simplex, log = TRUE)

    denominator <- target_log_dist(current_tau, current_simplex) +
        ddirichlet(proposed_simplex, current_simplex, log = TRUE)

    log_ratio <- numerator - denominator
    if (log_ratio > log(runif(1))) {
        current_simplex <- proposed_simplex
        acceptance_simplex[m] <- 1
    }
    simplexes[m, ] <- current_simplex
    ####
    #scale_tau <- update_scale(scale_tau, mean(acceptance_tau), target_acc = 0.55, it = 3)
    #scale_simplex <- update_scale(scale_simplex, mean(acceptance_simplex),
    #                              target_acc = 0.33, it = 3, c1 = 0.8, c0 = -1)
    #print(c(scale_tau = scale_tau, scale_simplex = scale_simplex))
}

mean(acceptance_tau[-seq_len(500L)])
mean(acceptance_simplex[-seq_len(500L)])

####

plot(taus[-seq_len(500L)], type = "l")

simplexes <- simplexes[-seq_len(500L), ]
plot(simplexes[, 1], type = "l")
plot(simplexes[, 2], type = "l")
plot(simplexes[, 3], type = "l")
plot(simplexes[, 4], type = "l")
plot(simplexes[, 5], type = "l")
plot(simplexes[, 6], type = "l")

####

mean_tau <- mean(taus[-seq_len(500L)])
mean_simplex <- colMeans(simplexes)

cor2cov(R, p * mean_tau * mean_simplex)
D

##########################################################################################
##########################################################################################

# sampling using the half-t approach

p <- ncol(D)
R <- cov2cor(D)
sds <- sqrt(diag(D))

D <- cor2cov(R, sds = sds)

b <- MASS::mvrnorm(1000, rep(0, p), D)

target_log_dist <- function (sds) {
    p <- length(sds)
    D <- cor2cov(R, sds^2)
    log_p_b <- sum(dmvnorm(b, rep(0, p), D, log = TRUE, prop = FALSE))
    log_p_tau <- sum(dht(sds, sigma = 15, df = 3, log = TRUE))
    log_p_b + log_p_tau
}

M <- 4000
acceptance_sds <- res_sds <- matrix(0.0, M, p)
current_sds <- sds
scale_sds <- rep(0.05, p)
if (p > 4)
    scale_sds[3:4] <- c(0.011)
#scale_sds <- 0.1 / sds
for (m in seq_len(M)) {
    for (i in seq_len(p)) {
        current_sds_i <- current_sds[i]
        scale_sds_i <- scale_sds[i]
        log_mu_i <- log(current_sds_i) - 0.5 * scale_sds_i^2
        proposed_sds_i <- rlnorm(1L, log_mu_i, scale_sds_i)
        pr <- current_sds
        pr[i] <- proposed_sds_i
        numerator_i <- target_log_dist(pr) +
            dlnorm(current_sds_i, log_mu_i, scale_sds_i, log = TRUE)
        denominator_i <- target_log_dist(current_sds) +
            dlnorm(proposed_sds_i, log_mu_i, scale_sds_i, log = TRUE)
        log_ratio_i <- numerator_i - denominator_i
        if (log_ratio_i > log(runif(1))) {
            current_sds <- pr
            acceptance_sds[m, i] <- 1
        }
        res_sds[m, i] <- current_sds[i]
    }
}

colMeans(acceptance_sds[-seq_len(1000L), ])

####

res_sds <- res_sds[-seq_len(1000L), ]
plot(res_sds[, 1], type = "l")
plot(res_sds[, 2], type = "l")
plot(res_sds[, 3], type = "l")
plot(res_sds[, 4], type = "l")
plot(res_sds[, 5], type = "l")
plot(res_sds[, 6], type = "l")

####

mean_sds <- colMeans(res_sds)

cor2cov(R, sds = mean_sds)
D





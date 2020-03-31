library("survival")
library("nlme")
library("GLMMadaptive")
library("splines")
#library("Formula")
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

cor2cov <- function (R, vars) {
    p <- nrow(R)
    if (length(vars) != p) {
        stop("incorrect length of 'vars' argument.")
    }
    sds <- sqrt(vars)
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
    #exp(log(x) - log(rowSums(x)))
}

##########################################################################################
##########################################################################################

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
    log_p_b <- sum(dmvnorm(b, rep(0, p), D, log = TRUE))
    log_p_tau <- dgamma(tau, 1, 1, log = TRUE)
    log_p_simplex <- ddirichlet(simplex, rep(1, p), log = TRUE)
    log_p_b + log_p_tau + log_p_simplex
}

M <- 2000
acceptance_tau <- taus <- numeric(M)
acceptance_simplex <- numeric(M)
simplexes <- matrix(0.0, M, p)
current_tau <- tau
current_simplex <- simplex
scale_tau <- 0.03
scale_simplex <- 1e7
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
    #if (FALSE) {
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
    #}
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







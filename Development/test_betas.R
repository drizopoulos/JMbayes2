library("nlme")
data("pbc2", package = "JMbayes2")
data("pbc2.id", package = "JMbayes2")
fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2,
           random = ~ year | id)

# response vector and design matrices
id <- pbc2$id
id <- match(id, unique(id))
n <- length(unique(id))
y <- log(pbc2$serBilir)
X <- model.matrix(~ year * sex, data = pbc2)
Z <- model.matrix(~ year, data = pbc2)
X_dot <- rbind(c(1, 0, 1, 0), c(0, 1, 0, 1))
X_dot <- X_dot[rep(1:2, n), ]
X_dot[, 3:4] <- X_dot[, 3:4] * rep(pbc2.id$sex == "female", each = 2L)
subjs_ind <- rep(seq_len(n), each = 2L)

# initial values
betas <- fixef(fm1)
betas[] <- 0
b <- data.matrix(ranef(fm1))
D <- getVarCov(fm1)
invD <- solve(D)
sigma <- fm1$sigma

# priors
mean_betas <- rep(0.0, length(betas))
Tau_betas <- diag(0.01, length(betas))
Tau_mean_betas <- c(Tau_betas %*% mean_betas)

# scale factor for proposal of random effects
scales <- 0.1 + 0 * b

robbins_monro <- function (scale, acceptance, it, target_acceptance = 0.45) {
    step_length <- scale / (target_acceptance * (1 - target_acceptance))
    ifelse(acceptance,
           scale + step_length * (1 - target_acceptance) / it,
           scale - step_length * target_acceptance / it)
}

# We will only update betas and b; for the betas we will use HC,
# but not for the b (i.e., the same as we do in JMbaye2)
M <- 3000L
res_betas <- matrix(0.0, M, length(betas))
for (m in seq_len(M)) {
    # update betas using Gibbs
    mean_u <- c(X_dot %*% betas)
    u <- c(t(b)) + mean_u
    sumXDu <- rep(0.0, 4)
    sumXDX <- matrix(0.0, 4, 4)
    for (i in seq_len(n)) {
        X_dot_i <- X_dot[subjs_ind == i, ]
        cprod <- crossprod(X_dot_i, invD)
        sumXDX <- sumXDX + cprod %*% X_dot_i
        sumXDu <- sumXDu + c(cprod %*% u[subjs_ind == i])
    }
    Sigma1 <- solve(sumXDX + Tau_betas)
    mean1 <- Sigma1 %*% (sumXDu + Tau_mean_betas)
    res_betas[m, ] <- betas <- MASS::mvrnorm(mu = mean1, Sigma = Sigma1)
    mean_u <- c(X_dot %*% betas)
    b <- matrix(u - mean_u, n, 2L, byrow = TRUE)
    # update b using Metropolis-Hastings
    mu <- c(X %*% betas) + rowSums(Z * b[id, ])
    log_denom <- c(rowsum(dnorm(y, mu, sigma, log = TRUE), id)) +
        JMbayes:::dmvnorm(b, rep(0, 2), invSigma = invD, log = TRUE)
    # loop over the random effects, i.e., update intercepts and slopes
    # separately, as we do in JMbayes2
    for (j in 1:2) {
        b_proposed <- b
        b_proposed[, j] <- rnorm(n, b[, j], scales[, j])
        mu_proposed <- c(X %*% betas) + rowSums(Z * b_proposed[id, ])
        log_numer <- c(rowsum(dnorm(y, mu_proposed, sigma, log = TRUE), id)) +
            JMbayes:::dmvnorm(b_proposed, rep(0, 2), invSigma = invD, log = TRUE)
        ratio <- exp(log_numer - log_denom)
        u <- runif(n)
        b[ratio > u, j] <- b_proposed[ratio > u, j]
        mu <- c(X %*% betas) + rowSums(Z * b[id, ])
        log_denom <- c(rowsum(dnorm(y, mu, sigma, log = TRUE), id)) +
            JMbayes:::dmvnorm(b, rep(0, 2), invSigma = invD, log = TRUE)
        # update scales using Robbins-Monro
        if (m > 20) {
            scales[, j] <- robbins_monro(scales[, j], ratio > u, m)
        }
    }
}

res_betas <- res_betas[-(1:1000), ]
cbind(Mean = colMeans(res_betas),
      t(apply(res_betas, 2L, quantile, probs = c(0.025, 0.975))))[, c(2, 1, 3)]
intervals(fm1)$fixed

plot(res_betas[, 1], type = "l")
plot(res_betas[, 2], type = "l")
plot(res_betas[, 3], type = "l")
plot(res_betas[, 4], type = "l")





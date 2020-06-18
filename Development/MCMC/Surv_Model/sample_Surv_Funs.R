log_density_surv <- function (bs_gammas, gammas, alphas) {
    lambda_H <- W0_H %*% bs_gammas + W_H %*% gammas
    for (i in seq_along(Wlong_H)) {
        lambda_H <- lambda_H + Wlong_H[[i]] %*% alphas[[i]]
    }
    lambda_h <- matrix(0.0, n, 1)
    if (length(which_event)) {
        lambda_h <- W0_h %*% bs_gammas + W_h %*% gammas
        for (i in seq_along(Wlong_h)) {
            W_h_i <- Wlong_h[[i]]
            lambda_h <- lambda_h + W_h_i %*% alphas[[i]]
        }
    }
    lambda_H2 <- matrix(0.0, nrow(Wlong_H2[[1]]), 1)
    if (length(which_interval)) {
        lambda_H2 <- W0_H2 %*% bs_gammas + W_H2 %*% gammas
        for (i in seq_along(Wlong_H2)) {
            W_H2_i <- Wlong_H2[[i]]
            lambda_H2 <- lambda_H2 + W_H2_i %*% alphas[[i]]
        }
    }
    H <- rowsum(exp(log_Pwk + lambda_H), group = id_H[[1]], reorder = FALSE)
    log_Lik_surv <- numeric(n)
    which_right_event <- c(which_right, which_event)
    if (length(which_right_event)) {
        log_Lik_surv[which_right_event] <- - H[which_right_event]
    }
    if (length(which_event)) {
        log_Lik_surv[which_event] <- log_Lik_surv[which_event] + lambda_h[which_event]
    }
    if (length(which_left)) {
        log_Lik_surv[which_left] <- log1p(- exp(- H[which_left]))
    }
    if (length(which_interval)) {
        H2 <- rowsum(exp(log_Pwk2 + lambda_H2), group = id_H2[[1]], reorder = FALSE)
        log_Lik_surv[which_interval] <- log(exp(- H[which_interval]) -
                                                exp(- (H2[which_interval] + H[which_interval])))
    }
    sum(log_Lik_surv, na.rm = TRUE)
}

rmvnorm <- function (n, mu = NULL, Sigma) {
    if (is.list(Sigma)) {
        ev <- Sigma$values
        evec <- Sigma$vectors
    } else {
        ed <- eigen(Sigma, symmetric = TRUE)
        ev <- ed$values
        evec <- ed$vectors
    }
    p <- length(ev)
    X <- tcrossprod(evec * rep(sqrt(pmax(ev, 0)), each = p),
                    matrix(rnorm(n * p), n))
    if (!is.null(mu))
        X <- drop(mu) + X
    X <- if (n == 1L) drop(X) else t.default(X)
    X
}

logPrior <- function (theta, mean_theta, Tau_theta, tau_theta) {
    z <- theta - mean_theta
    -0.5 * tau_theta * c(crossprod(z, Tau_theta) %*% z)
}

logPC_surv <- function (bs_gammas, gammas, alphas, tau_bs_gammas) {
    log_density_surv(bs_gammas, gammas, alphas) +
        logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas,
                 tau_bs_gammas)
}

robbins_monro_univ <- function (scale, acceptance_it, it, target_acceptance = 0.45) {
    step_length <- scale / (target_acceptance * (1 - target_acceptance))
    if (acceptance_it) {
        scale + step_length * (1 - target_acceptance) / it
    } else {
        scale - step_length * target_acceptance / it
    }
}


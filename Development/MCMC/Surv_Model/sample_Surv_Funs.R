log_density_surv2 <- function (W0H_bs_gammas, WH_gammas, WlongH_alphas,
                               W0h_bs_gammas, Wh_gammas, Wlongh_alphas,
                               W0H2_bs_gammas, WH2_gammas, WlongH2_alphas) {
    lambda_H <- W0H_bs_gammas + WH_gammas + WlongH_alphas
    lambda_h <- matrix(0.0, n, 1)
    if (length(which_event)) {
        lambda_h <- W0h_bs_gammas + Wh_gammas + Wlongh_alphas
    }
    lambda_H2 <- matrix(0.0, nrow(Wlong_H2[[1]]), 1)
    if (length(which_interval)) {
        lambda_H2 <- W0H2_bs_gammas + WH2_gammas + WlongH2_alphas
    }
    H <- group_sum(exp(log_Pwk + lambda_H), indFast_H)
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
        H2 <- group_sum(exp(log_Pwk2 + lambda_H2), indFast_H2)
        log_Lik_surv[which_interval] <- - H[which_interval] +
            log(- expm1(- H2[which_interval]))
    }
    sum(log_Lik_surv, na.rm = TRUE)
}

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

logPC_surv2 <- function (bs_gammas, gammas, alphas, tau_bs_gammas,
                         W0H_bs_gammas, WH_gammas, WlongH_alphas,
                         W0h_bs_gammas, Wh_gammas, Wlongh_alphas,
                         W0H2_bs_gammas, WH2_gammas, WlongH2_alphas) {
    log_density_surv2(W0H_bs_gammas, WH_gammas, WlongH_alphas,
                      W0h_bs_gammas, Wh_gammas, Wlongh_alphas,
                      W0H2_bs_gammas, WH2_gammas, WlongH2_alphas) +
        logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas,
                 tau_bs_gammas)
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

robbins_monro_mv2 <- function (scale, acceptance_it, it, dim,
                              target_acceptance = 0.25) {
    A <- 1 - 1 / dim
    alpha <- - qnorm(target_acceptance  / 2)
    B <- 0.5 * sqrt(2 * pi) * exp(alpha^2 / 2) / alpha
    C <- 1 / (dim * target_acceptance * (1 - target_acceptance))
    step_length <- scale * (A * B + C)
    den <- if (it > 299) max(299, it / dim) else it
    if (acceptance_it) {
        scale + step_length * (1 - target_acceptance) / den
    } else {
        scale - step_length * target_acceptance / den
    }
}

Wlong_alphas_fun <- function (Wlong, alphas) {
    out <- numeric(nrow(Wlong[[1L]]))
    for (i in seq_along(Wlong)) {
        out <- out + Wlong[[i]] %*% alphas[[i]]
    }
    out
}

group_sum <- function (x, ind) {
    xx <- c(0, cumsum(x)[ind])
    xx[-1L] - xx[-length(xx)]
}

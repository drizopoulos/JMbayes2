dmvnorm_chol <- function (x, chol_Sigma = NULL, chol_inv_Sigma = NULL, log = FALSE) {
    # this is a stripped down version of the multivariate normal density with mean zero
    # based on the Cholesky factor of the variance-covariance matrix Sigma
    if (!is.matrix(x))
        x <- rbind(x)
    p <- ncol(x)
    if (is.null(chol_Sigma) && is.null(chol_inv_Sigma))
        stop("'chol_Sigma' or 'chol_inv_Sigma' must be given.")
    invSigma <- if (!is.null(chol_Sigma)) {
        tcrossprod(backsolve(r = chol_Sigma, x = diag(p)))
    } else {
        crossprod(chol_inv_Sigma)
    }
    logdetSigma <- if (!is.null(chol_Sigma)) {
        2 * determinant(chol_Sigma)$modulus
    } else {
        - 2 * determinant(chol_inv_Sigma)$modulus
    }
    quad <- 0.5 * rowSums((x %*% invSigma) * x)
    fact <- -0.5 * (p * log(2 * pi) + logdetSigma)
    if (log) as.vector(fact - quad) else as.vector(exp(fact - quad))
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

cor2cov <- function (R, vars, sds = NULL) {
    p <- nrow(R)
    if (is.null(sds)) sds <- sqrt(vars)
    sds * R * rep(sds, each = p)
}

robbins_monro_univ <- function (scale, acceptance_it, it, target_acceptance = 0.45) {
    step_length <- scale / (target_acceptance * (1 - target_acceptance))
    if (acceptance_it) {
        scale + step_length * (1 - target_acceptance) / it
    } else {
        scale - step_length * target_acceptance / it
    }
}

logPC_D_sds <- function (sds, t_inv_L, half_t_df = 3, half_t_mean) {
    # log posterior conditional for sds
    log_p_b <- sum(dmvnorm_chol(b, chol_inv_Sigma = t_inv_L * rep(1 / sds, each = p),
                                log = TRUE))
    # prior is a half Student's-t
    log_p_sds <- sum(dht(sds, sigma = half_t_mean, df = half_t_df,
                         log = TRUE))
    log_p_b + log_p_sds
}

logPC_D_L <- function (L, sds, eta_LKJ = 2) {
    p <- length(sds)
    # log posterior conditional for the L matrix
    log_p_b <- sum(dmvnorm_chol(b, chol_Sigma = L * rep(sds, each = p),
                                log = TRUE))
    # LKJ prior on the Cholesky factor; check the following link for more info:
    # https://mc-stan.org/docs/2_18/functions-reference/cholesky-lkj-correlation-distribution.html
    log_p_L <- sum((p - 2:p + 2 * eta_LKJ - 2) * log(L[diags2]))
    log_p_b + log_p_L
}

reconstr_D <- function (L, sds) {
    p <- length(sds)
    LL <- matrix(0.0, p, p)
    LL[upper.tri(LL, TRUE)] <- L
    cor2cov(crossprod(LL), sds = sds)
}

deriv_L <- function (L, i, sds, log_target, eps = 1e-06,
                     direction = c("backward", "forward")) {
    direction <- match.arg(direction)
    if (direction == "forward") {
        L_eps1 <- L
        ##
        L_eps1[upper_tri_ind][i] <- L_eps1[upper_tri_ind][i] * (1 + eps)
        ll1 <- L_eps1[seq(1, colmn_ind[i] - 1), colmn_ind[i]]
        ss1 <- sum(ll1 * ll1)
        L_eps1[colmn_ind[i], colmn_ind[i]] <- sqrt(1 - ss1)
        ##
        (logPC_D_L(L_eps1, sds) - log_target) / eps
    } else {
        L_eps2 <- L
        ##
        L_eps2[upper_tri_ind][i] <- L_eps2[upper_tri_ind][i] * (1 - eps)
        ll2 <- L_eps2[seq(1, colmn_ind[i] - 1), colmn_ind[i]]
        ss2 <- sum(ll2 * ll2)
        if (ss2 > 1) return(as.numeric(NA))
        L_eps2[colmn_ind[i], colmn_ind[i]] <- sqrt(1 - ss2)
        ##
        (log_target - logPC_D_L(L_eps2, sds)) / eps
    }
}

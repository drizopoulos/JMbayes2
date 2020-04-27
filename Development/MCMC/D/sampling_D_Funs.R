dmvnorm_chol <- function (x, mu, chol_Sigma = NULL, chol_inv_Sigma = NULL, log = FALSE) {
    if (!is.matrix(x))
        x <- rbind(x)
    if (is.matrix(mu)) {
        if (nrow(mu) != nrow(x))
            stop("incorrect dimensions for 'mu'.")
        p <- ncol(mu)
    } else {
        p <- length(mu)
        mu <- rep(mu, each = nrow(x))
    }
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
    ss <- x - mu
    quad <- 0.5 * rowSums((ss %*% invSigma) * ss)
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
    log_p_b <- sum(dmvnorm_chol(b, rep(0, p),
                                chol_inv_Sigma = t_inv_L * rep(1 / sds, each = p),
                                log = TRUE))
    # prior is a half Student's-t
    log_p_sds <- sum(dht(sds, sigma = half_t_mean, df = half_t_df,
                         log = TRUE))
    log_p_b + log_p_sds
}

logPC_D_L <- function (L, sds, eta_LKJ = 2) {
    # log posterior conditional for the L matrix
    diags <- L[diags]
    if (any(is.na(diags))) return(-Inf)
    test_PD_1 <- all(diags > 0)
    test_PD_2 <- all(abs(sqrt(colSums(L^2)) - 1) < sqrt(.Machine$double.eps))
    if (!test_PD_1 || !test_PD_2) return(-Inf)
    log_p_b <- sum(dmvnorm_chol(b, rep(0, p), chol_Sigma = L * rep(sds, each = p),
                                log = TRUE))
    # LKJ prior on the Cholesky factor
    log_p_L <- sum((p - 2:p + 2 * eta_LKJ - 2) * log(diags[-1L]))
    log_p_b + log_p_L
}




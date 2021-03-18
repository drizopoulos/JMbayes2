i_row <- function (x, i) x[i, ]

splt_REs <- function (b_mat, ind_RE) {
    n <- length(ind_RE)
    out <- vector("list", n)
    for (i in seq_len(n)) out[[i]] <- b_mat[, ind_RE[[i]], drop = FALSE]
    out
}

linpred_long <- function (X, betas, Z, b, id, level = 0) {
    out <- vector("list", length(X))
    for (i in seq_along(X)) {
        X_i <- X[[i]]
        Z_i <- Z[[i]]
        betas_i <- betas[[i]]
        b_i <- b[[i]]
        id_i <- id[[i]]
        out[[i]] <- c(X_i %*% betas_i)
        if (level > 0) {
            out[[i]] <- out[[i]] + as.vector(rowSums(Z_i * b_i[id_i, ]))
        }
    }
    out
}

mu_fun <- function (eta, link) {
    switch (link, "identity" = eta, "inverse" = 1 / eta, "logit" = plogis(eta),
            "probit" = pnorm(eta), "cloglog" = - exp(- exp(eta)) + 1.0,
            "log" = exp(eta))
}

rowQuantile <- function (m, probs, na.rm = TRUE) {
    apply(m, 1L, quantile, probs = probs, na.rm = na.rm)
}



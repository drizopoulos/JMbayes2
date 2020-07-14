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

rmvnorm_array <- function (n, S, sigmas) {
  out <- array(0.0, dim = c(n, dim(S)[2], dim(S)[3]))
  for(i in 1:dim(S)[3]) {
    A <- chol(S[, ,i])
    z <- matrix(rnorm(n * ncol(S[, ,i])), nrow = n, ncol = ncol(S[, ,i]))
    out[, ,i] <- sqrt(sigmas[i]) * z %*% t(A)
  }
  out
}

robbins_monro_univ <- function (scale, acceptance_it, it, target_acceptance = 0.45) {
  step_length <- scale / (target_acceptance * (1 - target_acceptance))
  if (acceptance_it) {
    scale + step_length * (1 - target_acceptance) / it
  } else {
    scale - step_length * target_acceptance / it
  }
}

# this can probably be coded more efficiently but not needed since
# this function can be removed when translated to c++ in the final code
to_list_b <- function(b, b.cols) {
  out <- vector('list', length = length(b.cols))
  cumsum_b.cols <- cumsum(b.cols)
  if (b.cols[1] > 1) {
    out[[1]] <- t(b[, 1:b.cols[1], ])
  } else {
    out[[1]] <- matrix(b[, b.cols[1], ], ncol = 1)
  }
  for (i in 2:length(b.cols)) {
    if (b.cols[i] > 1) {
      out[[i]] <- t(b[, (b.cols[i-1]+1):(cumsum_b.cols[i]), ])
    } else {
      out[[i]] <- matrix(b[, b.cols[i], ], ncol = 1)
    }
  }
  out
}

log_post_b <- function(X, betas, Z, b, b.cols, id, 
                       y, log_sigmas, Funs, mu_funs, nY, unq_idL, idL, 
                       D, 
                       bs_gammas, gammas, 
                       alphas, 
                       n, bnew) {
  #b_lst <- list(t(b[, 1:2, ]), t(b[, 3:4, ]), matrix(b[, 5, ], ncol = 1), matrix(b[, 6, ], ncol = 1))
  #b_lst <- list(t(b[, 1:2, ]), matrix(b[, 3, ], ncol = 1), matrix(b[, 4, ], ncol = 1))
  b_lst <- to_list_b(b, b.cols)
  linear_predictor <- linpred_mixed(X, betas, Z, b_lst, id)
  log_pyb <- log_density_mixed(y, linear_predictor, log_sigmas, Funs, mu_funs, nY, unq_idL, idL)
  log_pb <- dmvnorm(t(b[1, , ]), mu = rep(0, ncol(b)), Sigma = D, log = TRUE, prop = TRUE)
  log_ptb <- log_density_surv(bs_gammas = bs_gammas, gammas = gammas, alphas = alphas)
  log_pyb + log_pb + log_ptb
}

log_post_b_HC <- function(X, betas, Z, b, b.cols, id, 
                          y, log_sigmas, Funs, mu_funs, nY, unq_idL, idL, 
                          D, Xhc, columns_HC,
                          bs_gammas, gammas, 
                          alphas, 
                          n, bnew) {
  #b_lst <- list(t(b[, 1:2, ]), t(b[, 3:4, ]), matrix(b[, 5, ], ncol = 1), matrix(b[, 6, ], ncol = 1))
  #b_lst <- list(t(b[, 1:2, ]), matrix(b[, 3, ], ncol = 1), matrix(b[, 4, ], ncol = 1))
  b_lst <- to_list_b(b, b.cols)
  linear_predictor <- linpred_mixed(X, betas, Z, b_lst, id)
  log_pyb <- log_density_mixed(y, linear_predictor, log_sigmas, Funs, mu_funs, nY, unq_idL, idL)
  mean_b_i <- calculate_mean_b_i(Xhc, columns_HC, betas, b_lst, unq_idL)
  mean_b_i <- do.call(cbind, mean_b_i)
  log_pb <- dmvnorm(t(b[1, , ]), mu = mean_b_i, Sigma = D, log = TRUE, prop = TRUE)
  log_ptb <- log_density_surv(bs_gammas = bs_gammas, gammas = gammas, alphas = alphas)
  log_pyb + log_pb + log_ptb
}

create_HC_X <- function (TermsX, TermsZ, x, z, id, mfHC) {
  # function that creates the hierarchical centering version of the
  # design matrix for the fixed effects
  find_positions <- function (nams1, nams2) {
    nams1 <- gsub("^", "\\^", nams1, fixed = TRUE)
    vals <- c(glob2rx(nams1), glob2rx(paste0(nams1, ":*")),
              glob2rx(paste0("*:", nams1)))
    out <- sort(unique(unlist(lapply(vals, grep, x = nams2))))
    out
  }
  check_td <- function (x, id) {
    !all(sapply(split(x, id), function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
  }
  grep2 <- function (x, nams) grep(x, nams, fixed = TRUE)
  terms.labs_X <- attr(TermsX, "term.labels")
  terms.labs_Z <- attr(TermsZ, "term.labels")
  # check for time-varying covariates
  timeTerms <- if (length(terms.labs_Z)) {
    unlist(lapply(terms.labs_Z, grep2, nams = colnames(x)))
  }
  which_td <- unname(which(apply(x, 2, check_td, id = id)))
  all_TDterms <- unique(c(timeTerms, which_td))
  baseline <- seq_len(ncol(x))[-all_TDterms]
  ind_colmns <- c(list(baseline), lapply(colnames(z)[-1L], find_positions,
                                         nams2 = colnames(x)))
  ind_colmns2 <- seq_len(ncol(x))
  ind_colmns2 <- ind_colmns2[!ind_colmns2 %in% unlist(ind_colmns)]
  mfHC <- mfHC[!duplicated(id), ]
  Xhc <- if (length(terms.labs_Z)) {
    which.timevar <- unique(unlist(lapply(terms.labs_Z, grep2, nams = names(mfHC))))
    mfHC[which.timevar] <- lapply(mfHC[which.timevar],
                                  function (x) { x[] <- 1; x })
    model.matrix(TermsX, mfHC)
  } else {
    model.matrix(TermsX, mfHC)
  }
  index <- numeric(ncol(Xhc))
  for (i in seq_along(ind_colmns)) {
    index[ind_colmns[[i]]] <- i
  }
  list(Xhc = Xhc, columns_HC = index, columns_nHC = ind_colmns2)
}

calculate_u <- function (Xhc, columns_HC, betas, b, unq_idL) {
  u <- b
  for (i in seq_along(Xhc)) {
    Xhc_i <- Xhc[[i]]
    columns_HC_i <- columns_HC[[i]]
    betas_i <- betas[[i]]
    b_i <- b[[i]]
    unq_idL_i <- unq_idL[[i]]
    mean_b_i <- b_i * 0
    for (j in seq_len(ncol(b_i))) {
      index <- columns_HC_i == j
      mean_b_i[unq_idL_i, j] <- c(Xhc_i[, index, drop = FALSE] %*% betas_i[index])
    }
    u[[i]] <- b_i + mean_b_i
  }
  u
}

calculate_mean_b_i <- function (Xhc, columns_HC, betas, b, unq_idL) {
  u <- b
  for (i in seq_along(Xhc)) {
    Xhc_i <- Xhc[[i]]
    columns_HC_i <- columns_HC[[i]]
    betas_i <- betas[[i]]
    b_i <- b[[i]]
    unq_idL_i <- unq_idL[[i]]
    mean_b_i <- b_i * 0
    for (j in seq_len(ncol(b_i))) {
      index <- columns_HC_i == j
      mean_b_i[unq_idL_i, j] <- c(Xhc_i[, index, drop = FALSE] %*% betas_i[index])
    }
    u[[i]] <- mean_b_i
  }
  u
}



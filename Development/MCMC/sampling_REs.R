# Source Functions
source(file.path(getwd(), "R/jm.R"))
source(file.path(getwd(), "R/help_functions.R"))
source(file.path(getwd(), "Development/jm/R_to_Cpp.R"))
source(file.path(getwd(), "Development/jm/PBC_data.R"))

# load data
load(file = file.path(getwd(), "/Dev_Local/sample_case_env_02042020.RData"))

# note that for missing subjects for a long outcome the vector of random-effects
# equals zero
# lapply(idL, FUN = function(x) length(unique(x)))
# lapply(b, nrow)
# which(!unique(idL[[1]]) %in% unique(idL[[2]]))
# which(b[[2]][, 2] == 0)

#---------------------------------------------------------------
#                        FUNCTIONS
#---------------------------------------------------------------
# compund assignment operator += similar to c++
`%+=%` <- function(x1, x2) eval.parent(substitute(x1 <- x1 + x2))

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

# sample random values from multivariate normal distribution
mvrnorm_gp <- function(n, S, sigma) {
  A <- chol(S)
  z <- matrix(rnorm(n * ncol(S)), nrow = n, ncol = ncol(S))
  Y <- sqrt(sigma) * z %*% t(A)
  Y
}

mvrnorm_gp_array <- function (n, S, sigmas) {
  out <- array(0.0, dim = c(n, dim(S)[2], dim(S)[3]))
  for(i in 1:dim(S)[3]) {
    A <- chol(S[, ,i])
    z <- matrix(rnorm(n * ncol(S[, ,i])), nrow = n, ncol = ncol(S[, ,i]))
    out[, ,i] <- sqrt(sigmas[i]) %*% z %*% t(A)
  }
}

#log_long_density <- function(y, eta, families, sigmas, id, n) {
#  n_outcomes <- length(y)
#  out <- numeric(length = n)
#  for (i in 1:n_outcomes) {
#    id_i <- id[[i]]
#    y_i <- y[[i]]
#    eta_i <- eta[[i]]
#    if (families[[i]]$family == "gaussian") {
#      sigma_i <- sigmas[[i]]
#      log_dens = -0.5 * ((y_i - eta_i) / sigma_i)^2
#      out %+=% rowsum(log_dens, id_i)
#    } else if (families[[i]]$family == "binomial") {
#      if (families[[i]]$link == 'logit') {
#        pr = exp(eta_i) / (1 + exp(eta_i))
#        log_dens = y_i * log(pr) + (1 - y_i) * log(1 - pr)
#        out %+=% rowsum(log_dens, id_i)
#      }
#    }
#  }
#  out
#}


target_log_dist <- function(X, betas, Z, b, id, 
                            y, log_sigmas, Funs, mu_funs, nY, unq_idL, idL, 
                            D,  
                            log_Lik_surv) {
  #b_mat <- do.call(cbind, b)
  b_mat <- b
  b <- list(b[, 1:2], b[, 3:4], matrix(b[, 5], ncol = 1), matrix(b[, 6], ncol = 1))
  linear_predictor <- linpred_mixed(X, betas, Z, b, idL)
  log_pyb <- sum(log_density_mixed(y, linear_predictor, log_sigmas, Funs, mu_funs, nY, unq_idL, idL))
  log_pb <- sum(dmvnorm(b_mat, mu = rep(0, ncol(b_mat)), Sigma = D, log = TRUE, prop = FALSE))
  #Wlong_h_mat <- do.call(cbind, Wlong_h)
  #Wlong_H_mat <- do.call(cbind, Wlong_H)
  #alphas_vec <- do.call(c, alphas)
  #log_h <- W0_h %*% bs_gammas + W_h %*% gammas + Wlong_h_mat * alphas_vec
  #H <- rowSums(P * exp(W0_H %*% bs_gammas + Wlong_H_mat %*% alphas_vec))
  #delta[delta > 1] <- 1
  #log_ptb <- sum(log_Lik_surv[!is.na(log_Lik_surv) & log_Lik_surv != -Inf])
  log_pyb + log_ptb + log_pb
  log_pyb + log_pb
}

target_log_dist_2 <- function(b, log_pyb) {
  log_pb <- sum(dmvnorm(b, mu = rep(0, ncol(b)), Sigma = D, log = TRUE, prop = FALSE))
  log_pyb + log_pb
}
  


# sample random values from t distribution

# MCMC
M <- 3000
b.rows <- max(do.call(c, lapply(b, nrow)))
b.cols <- do.call(c, lapply(b, ncol))
bs <- array(0.0, dim = c(length(unq_idL[[1]]), sum(b.cols), M))
current_b <- do.call(cbind, b)
acceptance_b <- numeric(M)

for (m in seq_len(M)) {
  proposed_b <- mvrnorm_gp(b.rows, D, 1e-10)
  numerator <- target_log_dist(X, betas, Z, proposed_b, idL, 
                               y, log_sigmas, Funs, mu_funs, nY, unq_idL, idL, 
                               D, log_Lik_surv)
  denominator <- target_log_dist(X, betas, Z, current_b, idL, 
                                 y, log_sigmas, Funs, mu_funs, nY, unq_idL, idL, 
                                 D, log_Lik_surv)
  log_ratio <- numerator - denominator
  if (log_ratio > log(runif(1))) {
    current_b <- proposed_b
    acceptance_b[m] <- 1
  }
  bs[, ,m] <- current_b
}

for (m in seq_len(M)) {
  proposed_b <- mvrnorm_gp(b.rows, D, 1e-3)
  numerator <- target_log_dist_2(proposed_b, log_pyb)
  denominator <- target_log_dist_2(current_b, log_pyb)
  log_ratio <- numerator - denominator
  if (log_ratio > log(runif(1))) {
    current_b <- proposed_b
    acceptance_b[m] <- 1
  }
  bs[, ,m] <- current_b
}


mean(acceptance_b[-seq_len(500L)])
plot(bs[2, 3, ], type = 'l')

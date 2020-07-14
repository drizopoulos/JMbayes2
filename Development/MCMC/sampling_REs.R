#########################################################################################################
# SIMULATE DATA FROM A JOINT MODEL - SAVE THE RANDOM-EFFECTS - AND CHECK WHETHER IT WORKS - 3TO5 DATASETS
# TO IMPLEMENT HIERARCHICAL CENTERING
# CENTERING OF THE X AND Z MATRICES
# CHECK SYSTEM TIME
# MULTIPLE CHAINS (JITTER (NOISE FROM UNIFORM) OR NORMAL DISTRIBUTION WITH VARIANCE)
#########################################################################################################

library("survival")
library("nlme")
library("GLMMadaptive")
library("splines")
data("pbc2", package = "JM")
data("pbc2.id", package = "JM")

# Source Functions
source(file.path(getwd(), "R/jm.R"))
source(file.path(getwd(), "R/help_functions.R"))
source(file.path(getwd(), "Development/jm/R_to_Cpp.R"))
source(file.path(getwd(), "Development/jm/PBC_data.R"))
load(file = file.path(getwd(), "/Development/Function_Inputs_R_to_Cpp/init_vals_surv.RData"))
load(file = file.path(getwd(), "/Development/Function_Inputs_R_to_Cpp/init_surv.RData"))

# load data
if (length(grep('gpapageorgiou', getwd())) == 1) {
  load(file = file.path(getwd(), "/Dev_Local/sample_case_env_02042020.RData"))
  load(file = file.path(getwd(), "/Dev_Local/sample_case_env_testjm_08042020.RData"))
} else {
  source(file.path(getwd(), "Development/jm/sample_case.R"))
  source(file.oath(getwd(), "Development/jm/test_jm.R"))
}


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
    out[, ,i] <- sqrt(sigmas[i]) * z %*% t(A)
  }
  out
}

log_dens_surv <- function (bs_gammas, n, Data, Wlong_h, Wlong_H, gammas, alphas, id_H) {
  W0_h <- Data$W0_h
  W0_H <- Data$W0_H
  #W0_H2 <- Data$W0_H2
  W_h <- Data$W_h
  W_H <- Data$W_H
  #W_H2 <- Data$W_H2
  log_Pwk <- Data$log_Pwk
  log_Pwk2 <- Data$log_Pwk2
  which_event <- Data$which_event
  #which_interval <- Data$which_interval
  which_left <- Data$which_left
  which_right <- Data$which_right
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
  #lambda_H2 <- matrix(0.0, nrow(Wlong_H2[[1]]), 1)
  #if (length(which_interval)) {
  #  lambda_H2 <- W0_H2 %*% bs_gammas + W_H2 %*% gammas
  #  for (i in seq_along(Wlong_H2)) {
  #    W_H2_i <- Wlong_H2[[i]]
  #    lambda_H2 <- lambda_H2 + W_H2_i %*% alphas[[i]]
  #  }
  #}
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
  #if (length(which_interval)) {
  #  H2 <- rowsum(exp(log_Pwk2 + lambda_H2), group = id_H2[[1]], reorder = FALSE)
  #  log_Lik_surv[which_interval] <- log(exp(- H[which_interval]) -
  #                                        exp(-H2[which_interval]))
  #}
  - sum(log_Lik_surv, na.rm = TRUE)
}

target_log_dist <- function(X, betas, Z, b, id, 
                            y, log_sigmas, Funs, mu_funs, nY, unq_idL, idL, 
                            D, 
                            Data, 
                            Wlong_h, Wlong_H,
                            bs_gammas, gammas, 
                            alphas, 
                            id_H, 
                            n, bnew) {
  b_lst <- list(t(b[, 1:2, ]), t(b[, 3:4, ]), matrix(b[, 5, ], ncol = 1), matrix(b[, 6, ], ncol = 1))
  linear_predictor <- linpred_mixed(X, betas, Z, b_lst, id)
  log_pyb <- log_density_mixed(y, linear_predictor, log_sigmas, Funs, mu_funs, nY, unq_idL, idL)
  log_pb <- dmvnorm(t(b[1, , ]), mu = rep(0, ncol(b)), Sigma = D, log = TRUE, prop = TRUE)
  #log_pb <- dmvnorm(t(b[1, , ]), mu = t(bnew[1, , ]), Sigma = D, log = TRUE, prop = FALSE)
  log_ptb <- log_dens_surv(bs_gammas = bs_gammas, n = n, Data = Data, Wlong_h = Wlong_h, Wlong_H = Wlong_H, 
                           gammas = gammas, alphas = alphas, id_H = id_H)
  log_pyb + log_pb + log_ptb
}

robbins_monro_univ <- function (scale, acceptance_it, it, target_acceptance = 0.45) {
  step_length <- scale / (target_acceptance * (1 - target_acceptance))
  if (acceptance_it) {
    scale + step_length * (1 - target_acceptance) / it
  } else {
    scale - step_length * target_acceptance / it
  }
}

#

unq_idL <- lapply(idL, unique)

# MCMC
M <- 10000
b.rows <- max(do.call(c, lapply(b, nrow)))
b.cols <- do.call(c, lapply(b, ncol))
bs <- array(0.0, dim = c(M, sum(b.cols), length(unq_idL[[1]])))
b_mat <- do.call(cbind, b)
current_b <- array(0.0, dim = c(1, sum(b.cols), length(unq_idL[[1]])))
for (i in 1:sum(b.cols)) {
  current_b[, i, ] <- b_mat[, i]
}
#init_b <- current_b
acceptance_b <- matrix(0.0, nrow = b.rows, ncol = M)
sigmas <- rep(6 / sum(b.cols), b.rows)
vcov_prop_RE <- test$vcov_prop$vcov_prop_RE
#proposed_b <- mvrnorm_gp_array(1, vcov_prop_RE, sigmas)
log_us_RE <- matrix(log(runif(b.rows * M)), nrow = b.rows, ncol = M)
#invD <- solve(D)


for (m in seq_len(M)) {
  proposed_b <- mvrnorm_gp_array(1, vcov_prop_RE, sigmas) + current_b
  numerator <- target_log_dist(X = X, betas = betas, Z = Z, b = proposed_b, id = idL_lp, 
                               y = y, log_sigmas = log_sigmas, Funs = Funs, mu_funs = mu_funs, nY = nY, unq_idL = unq_idL, idL = idL, 
                               D = D,
                               Data = Data,
                               Wlong_h = Wlong_h, Wlong_H = Wlong_H, 
                               bs_gammas = init_surv$bs_gammas, gammas = test$gammas, 
                               alphas = init_surv$alphas,
                               id_H = id_H, 
                               n = b.rows, bnew = proposed_b)
  denominator <- target_log_dist(X = X, betas = betas, Z = Z, b = current_b, id = idL_lp, 
                                 y = y, log_sigmas = log_sigmas, Funs = Funs, mu_funs = mu_funs, nY = nY, unq_idL = unq_idL, idL = idL, 
                                 D = D,
                                 Data = Data,
                                 Wlong_h = Wlong_h, Wlong_H = Wlong_H, 
                                 bs_gammas = init_surv$bs_gammas, gammas = test$gammas, 
                                 alphas = init_surv$alphas,
                                 id_H = id_H, 
                                 n = b.rows, bnew = current_b)
  log_ratio <- numerator - denominator
  for (i in 1:length(log_ratio)) {
    if (log_us_RE[i, m] < log_ratio[i]) {
      acceptance_b[i, m] <- 1.0
      current_b[, ,i] <- proposed_b[, ,i]
    }
    bs[m, ,i] <- current_b[, ,i]
    if (m > 20) {
      sigmas[i] <- robbins_monro_univ(scale = sigmas[i],
                                         acceptance_it = acceptance_b[i, m],
                                         it = m, target_acceptance = 0.4)
    }
  }
}


#mean(acceptance_b[312, ][-seq_len(500L)])
plot(bs[, 6, 1], type = 'l')
apply(acceptance_b, MARGIN = 1, mean)

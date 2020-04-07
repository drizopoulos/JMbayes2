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

# sample random values from multivariate normal distribution
mvrnorm_gp <- function(n, S) {
  A <- chol(S)
  z <- matrix(rnorm(n * ncol(S)), nrow = n, ncol = ncol(S))
  Y <- z %*% t(A)
  Y
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
                            Wlong_h, Wlong_H, 
                            alphas, 
                            W0_h, W0_H, W_h, delta, 
                            log_Lik_surv) {
  b_mat <- do.call(cbind, b)
  linear_predictor <- linpred_mixed(X, betas, Z, b, id)
  log_pyb <- log_density_mixed(y, linear_predictor, log_sigmas, Funs, mu_funs, nY, unq_idL, idL)
  log_pb <- JMbayes:::dmvnorm2(b_mat, mean = rep(0, ncol(b_mat)), sigma = D, logd = TRUE)
  #log_pb <- JMbayes:::dmvnorm(x = b_mat, mu = rep(0, ncol(b_mat)), Sigma = D, log = TRUE, prop = FALSE)
  #Wlong_h_mat <- do.call(cbind, Wlong_h)
  #Wlong_H_mat <- do.call(cbind, Wlong_H)
  #alphas_vec <- do.call(c, alphas)
  #log_h <- W0_h %*% bs_gammas + W_h %*% gammas + Wlong_h_mat * alphas_vec
  #H <- rowSums(P * exp(W0_H %*% bs_gammas + Wlong_H_mat %*% alphas_vec))
  #delta[delta > 1] <- 1
  log_ptb <- log_Lik_surv
  log_pyb + log_ptb + log_pb
}
  


# sample random values from t distribution

# MCMC
M <- 3000
b.rows <- max(do.call(c, lapply(b, nrow)))
b.cols <- do.call(c, lapply(b, ncol))
bs <- array(0.0, dim = c(length(idL[[1]]), sum(b.cols), M))


for (m in seq_len(M)) {
  accept_b <- matrix(0.0, nrow = b.rows, ncol = M)
  
  
}
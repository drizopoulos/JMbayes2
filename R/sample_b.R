#source(file.path(getwd(), "Development/MCMC/b/examples.R"))
source(file.path(getwd(), "Development/MCMC/b/examples_2.R"))
source(file.path(getwd(), "Development/MCMC/b/sample_b_FUNS.R"))

#----------------------------------
# WITHOUT HC
#----------------------------------
M <- 5000L
b.rows <- max(do.call(c, lapply(b, nrow)))
b.cols <- do.call(c, lapply(b, ncol))
bs <- array(0.0, dim = c(M, sum(b.cols), b.rows))
acceptance_b <- matrix(0.0, nrow = b.rows, ncol = M)
b_mat <- do.call(cbind, b)
current_b <- array(0.0, dim = c(1, sum(b.cols), length(unq_idL[[1]])))
for (i in 1:sum(b.cols)) {
  current_b[, i, ] <- b_mat[, i]
}
sigmas <- rep(0.6 / sum(b.cols), b.rows)
vcov_prop_RE <- test$vcov_prop$vcov_prop_RE
log_us_RE <- matrix(log(runif(b.rows * M)), nrow = b.rows, ncol = M)
bs_gammas <- test$initial_values$bs_gammas
gammas <- test$initial_values$gammas
alphas <- test$initial_values$alphas
D <- test$initial_values$D

system.time({
  for (m in seq_len(M)) {
    proposed_b <- rmvnorm_array(1, vcov_prop_RE, sigmas) + current_b
    numerator <- log_post_b(X = X, betas = betas, Z = Z, b = proposed_b, b.cols = b.cols, id = idL_lp, 
                            y = y, log_sigmas = log_sigmas, Funs = Funs, mu_funs = mu_funs, nY = nY, unq_idL = unq_idL, idL = idL, 
                            D = D,
                            bs_gammas = bs_gammas, gammas = gammas, 
                            alphas = alphas,
                            n = b.rows, bnew = proposed_b)
    denominator <- log_post_b(X = X, betas = betas, Z = Z, b = current_b, b.cols = b.cols, id = idL_lp, 
                              y = y, log_sigmas = log_sigmas, Funs = Funs, mu_funs = mu_funs, nY = nY, unq_idL = unq_idL, idL = idL, 
                              D = D,
                              bs_gammas = bs_gammas, gammas = gammas, 
                              alphas = alphas,
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
})

par(mfrow = c(4, 4))

for (i in 1:dim(bs)[2]) {
  for (j in 1:dim(bs)[3]) {
    plot(bs[, i, j], type = 'l')
  }
}

par(mfrow = c(1, 1))

plot(bs[, 1, 1], type = 'l')
apply(acceptance_b, MARGIN = 1, mean)



#--------------------
# WITH HC
#--------------------
M <- 5000L
b.rows <- max(do.call(c, lapply(b, nrow)))
b.cols <- do.call(c, lapply(b, ncol))
bs <- array(0.0, dim = c(M, sum(b.cols), b.rows))
acceptance_b <- matrix(0.0, nrow = b.rows, ncol = M)
u <- calculate_u(Xhc, columns_HC, betas, b, unq_idL)
u_mat <- do.call(cbind, u)
current_u <- array(0.0, dim = c(1, sum(b.cols), length(unq_idL[[1]])))
for (i in 1:sum(b.cols)) {
  current_u[, i, ] <- u_mat[, i]
}
sigmas <- rep(0.6 / sum(b.cols), b.rows)
vcov_prop_RE <- test$vcov_prop$vcov_prop_RE
log_us_RE <- matrix(log(runif(b.rows * M)), nrow = b.rows, ncol = M)
bs_gammas <- test$initial_values$bs_gammas
gammas <- test$initial_values$gammas
alphas <- test$initial_values$alphas
D <- test$initial_values$D

system.time({
  for (m in seq_len(M)) {
    proposed_u <- rmvnorm_array(1, vcov_prop_RE, sigmas) + current_u
    numerator <- log_post_b_HC(X = X, betas = betas, Z = Z, b = proposed_u, b.cols = b.cols, id = idL_lp, 
                            y = y, log_sigmas = log_sigmas, Funs = Funs, mu_funs = mu_funs, nY = nY, unq_idL = unq_idL, idL = idL, 
                            D = D, Xhc = Xhc, columns_HC = columns_HC,
                            bs_gammas = bs_gammas, gammas = gammas, 
                            alphas = alphas,
                            n = b.rows, bnew = proposed_u)
    denominator <- log_post_b_HC(X = X, betas = betas, Z = Z, b = current_u, b.cols = b.cols, id = idL_lp, 
                              y = y, log_sigmas = log_sigmas, Funs = Funs, mu_funs = mu_funs, nY = nY, unq_idL = unq_idL, idL = idL, 
                              D = D, Xhc = Xhc, columns_HC = columns_HC,
                              bs_gammas = bs_gammas, gammas = gammas, 
                              alphas = alphas,
                              n = b.rows, bnew = current_u)
    log_ratio <- numerator - denominator
    for (i in 1:length(log_ratio)) {
      if (log_us_RE[i, m] < log_ratio[i]) {
        acceptance_b[i, m] <- 1.0
        current_u[, ,i] <- proposed_u[, ,i]
      }
      bs[m, ,i] <- current_u[, ,i]
      if (m > 20) {
        sigmas[i] <- robbins_monro_univ(scale = sigmas[i],
                                        acceptance_it = acceptance_b[i, m],
                                        it = m, target_acceptance = 0.4)
      }
    }
  }
})

par(mfrow = c(4, 4))

for (i in 1:dim(bs)[2]) {
  for (j in 1:dim(bs)[3]) {
    plot(bs[, i, j], type = 'l')
  }
}

par(mfrow = c(1, 1))
plot(bs[, 1, 1], type = 'l')
apply(acceptance_b, MARGIN = 1, mean)





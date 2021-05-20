extract_last_iterations <- function(x) {
  x_mcmc <- x$mcmc
  n_chains <- x$control$n_chains
  nams_x <- names(x_mcmc)
  nams_special <- c(nams_x[grep('betas', nams_x)], 'b', 'D', 'sigmas')
  ind_RE <- x$model_data$ind_RE
  dim_D <- ncol(x$statistics$Mean$b)
  last_iter_x <- vector('list', n_chains)
  has_sigmas <- as.integer(x$initial_values$log_sigmas > -20)
  for (i in 1:n_chains) {
    last_iter_x[[i]] <- extract_mcmc_as_inits(x_mcmc, i = i, nams_x = nams_x, 
                                              nams_special = nams_special, ind_RE = ind_RE, dim_D, 
                                              has_sigmas = has_sigmas)
  }
  last_iter_x
}


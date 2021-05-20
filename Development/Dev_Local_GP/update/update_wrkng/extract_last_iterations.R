extract_last_iterations <- function(x) {
  x_mcmc <- x$mcmc
  n_chains <- x$control$n_chains
  namsX <- names(x_mcmc)
  nams_exclude <- c('betas', 'b', 'D')
  last_iter_b <- list(NULL)
  last_iter <- list(NULL)
  last_iter_D <- list(NULL)
  last_iter_betas <- list(NULL)
  for (i in 1:n_chains) {
    last_iter_b[[i]] <- as.list(lapply(x_mcmc[namsX %in% 'b'], function(x, i) x[[i]][, , dim(x[[i]])[length(dim(x[[i]]))], drop = FALSE], i = i))
    last_iter_D[[i]] <- lapply(x_mcmc[namsX %in% 'D'], function(x, i) x[[i]][nrow(x[[i]]), , drop = FALSE], i = i)
    last_iter_betas[[i]] <- as.list(lapply(x_mcmc[namsX %in% 'betas'], function(x, i) x[[i]][nrow(x[[i]]), , drop = TRUE], i = i))
    last_iter[[i]] <- lapply(x_mcmc[!namsX %in% nams_exclude], function(x, i) x[[i]][nrow(x[[i]]), , drop = TRUE], i = i)
    last_iter[[i]] <- c(last_iter_betas[[i]], last_iter[[i]], last_iter_D[[i]], last_iter_b[[i]])
  }
  last_iter
}


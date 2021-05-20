extract_mcmc_as_inits <- function(x_mcmc, i, nams_x, nams_special, ind_RE, dim_D, has_sigmas) {
  c(lapply(x_mcmc[!nams_x %in% nams_special], function(x, i) x[[i]][nrow(x[[i]]), , drop = TRUE], i = i), 
    'log_sigmas' = unname(lapply(x_mcmc[nams_x %in% 'sigmas'], function(x, i, has_sigmas) ifelse(has_sigmas, log(x[[i]][nrow(x[[i]]), , drop = TRUE]), -20), i = i, has_sigmas = has_sigmas)),
    list(betas = unname(lapply(x_mcmc[nams_x[grep('betas', nams_x)]], function(x, i) x[[i]][nrow(x[[i]]), , drop = TRUE], i = i))), 
    list(b = mapply(function(x, y, i) x[[i]][,y , dim(x[[i]])[length(dim(x[[i]]))], drop = TRUE], 
                    x_mcmc[nams_x %in% 'b'], ind_RE, 
                    i = i, SIMPLIFY = FALSE, USE.NAMES = FALSE)), 
    lapply(x_mcmc[nams_x %in% 'D'], function(x, i, dim_D) {
      D <- matrix(NA, ncol = dim_D, nrow = dim_D)
      D_vec <- x[[i]][nrow(x[[i]]), , drop = TRUE]
      D[upper.tri(D, diag = TRUE)] <- D_vec[order(names(D_vec))]
      D[lower.tri(D, diag = FALSE)] <- D[upper.tri(D, diag = FALSE)]
      D
    }, 
    i = i, dim_D = dim_D)
  )
}

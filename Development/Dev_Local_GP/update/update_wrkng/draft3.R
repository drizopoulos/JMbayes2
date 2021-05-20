library(parallel)

check1 <- function(chain, xlist) {
  xlist[[1]]
}

check2 <- function(chain, xlist) {
  xlist[[1]] <- 'changed'
  xlist
}

check3 <- function(chain, xlist) {
  out <- xlist[[1]] 
  out <- 'changed'
  out
}

all.equal(model_data_list[[1]], model_data[[2]])


xlist <- list('A' = 'A', 'B' = 'B', 'C' = 'C')
chains = 1:3

ncores <- detectCores() - 1
cores <- min(ncores, length(chains))
cl <- makeCluster(cores)
parallel::parLapply(cl, chains, check1, xlist)
parallel::parLapply(cl, chains, check2, xlist)
parallel::parLapply(cl, chains, check3, xlist)

stopCluster(cl)


if (is.null(initial_values)) {
  initial_values <- lapply(seq_len(con$n_chains), 
                           function(x, ...) list(betas = betas, log_sigmas = log_sigmas, D = D, 
                                                 b = b, bs_gammas = bs_gammas, gammas = gammas, 
                                                 alphas = alphas, tau_bs_gammas = rep(20, n_strata)), 
                           betas = betas, log_sigmas = log_sigmas, 
                           D = D, b = b, bs_gammas = bs_gammas, gammas = gammas, 
                           alphas = alphas, n_strata = n_strata)
}
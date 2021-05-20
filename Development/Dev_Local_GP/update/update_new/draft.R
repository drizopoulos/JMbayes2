initial_values

list(last_iterations[[1]]$betas1)

initial_values$
x_mcmc$sigmas$`1`
x_mcmc$b$`1`
x$mcmc$b$`1`

x <- joint_model_fit_2
?list
klain <- vector('list', 3)

extract_mcmc <- function(x_mcmc, i, nams_x, nams_special, ind_RE, dim_D) {
  
  c(lapply(x_mcmc[!nams_x %in% nams_special], function(x, i) x[[i]][nrow(x[[i]]), , drop = TRUE], i = i), 
    list(betas = unname(lapply(x_mcmc[nams_x[grep('betas', nams_x)]], function(x, i) x[[i]][nrow(x[[i]]), , drop = TRUE], i = i))), 
    list(b = mapply(function(x, y, i) x[[i]][,y , dim(x[[i]])[length(dim(x[[i]]))], drop = FALSE], 
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


chk <- x_mcmc
dim(x_mcmc$b$`1`)


JMbayes2:::mapply2(function(x, y, i) dim(x[[i]]), 
                   x_mcmc[names(x_mcmc) %in% 'b'], x$model_data$ind_RE, 
                   i = i)

?mapply
mapply(function(x, y, i) x[[i]][,y , dim(x[[i]])[length(dim(x[[i]]))], drop = FALSE], 
                   x_mcmc[names(x_mcmc) %in% 'b'], x$model_data$ind_RE, 
                   i = i, SIMPLIFY = FALSE, USE.NAMES = FALSE)

matrix(x_mcmc$D$`1`[1, ], nrow = ncol(x_mcmc$b$`1`), ncol = ncol(x_mcmc$b$`1`))

JMbayes2:::diag(x_mcmc$D$`1`[1, ])
initial_values$D

x$statistics$Mean$D

Matrix::Diagonal(4, x = x_mcmc$D$`1`[1, ])

JMbayes2:::.bdiag(x_mcmc$D$`1`)

upper.tri(x_mcmc$D$`1`[1, ])

?upper.tri
initial_values$D
?bdiag

last_iter_x[[1]]$D
x_mcmc$D$`1`[3000, ]

extract_last_iterations <- function(x) {
  x_mcmc <- x$mcmc
  n_chains <- x$control$n_chains
  nams_x <- names(x_mcmc)
  nams_special <- c(nams_x[grep('betas', nams_x)], 'b', 'D')
  ind_RE <- x$model_data$ind_RE
  dim_D <- ncol(x$statistics$Mean$b)
  last_iter_x <- vector('list', n_chains)
  
  for (i in 1:n_chains) {
    last_iter_x[[i]] <- extract_mcmc(x_mcmc, i = i, nams_x = nams_x, 
                                     nams_special = nams_special, ind_RE = ind_RE, dim_D)
  }
  last_iter_x
  
  
  
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



upper.tri(matrix(x_mcmc$D$`1`[1, ]), diag = TRUE)
demek <- matrix(NA, ncol = 4, nrow = 4)
demek[upper.tri(demek, diag = TRUE)] <- x_mcmc$D$`1`[1, ][order(names(x_mcmc$D$`1`[1, ]))]
demek[lower.tri(demek, diag = FALSE)] <- demek[upper.tri(demek, diag = FALSE)]

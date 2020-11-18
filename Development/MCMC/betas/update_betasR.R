update_betasR <- function(betas,
                          prior_mean_betas_HC, prior_Tau_betas_HC,
                          b_mat, D, X_dot,
                          ind_FE,
                          ind_FE_HC,
                          id_patt,
                          ind_RE_patt,
                          ind_FE_patt,
                          it,
                          y) {
  
  
  betas_vec <- unlist(betas, use.names = FALSE)
  
  n_b <- nrow(b_mat)
  q <- ncol(b_mat)
  patt_count <- length(unique(id_patt))
  p_HC <- length(ind_FE_HC) 
  sum_JXDXJ <- matrix(0, nrow= p_HC, ncol= p_HC)
  sum_JXDu <- rep(0, p_HC)
  D_inv <- list()
  
  for(i in seq_len(patt_count)) { # new 
      
    if(length(ind_RE_patt[[i]])==0) next
    
    D_inv[[i]] <- solve(D[ind_RE_patt[[i]], ind_RE_patt[[i]]])
    
  }
  
  
  for(i in seq_len(n_b)) {
    
    patt_i <- id_patt[i]
    ind_FE_i <- ind_FE_patt[[patt_i]]
    ind_RE_i <- ind_RE_patt[[patt_i]]
    
    X_dot_i <- X_dot[(1 + (i - 1) * q):(q + (i - 1) * q), ]
    X_dot_i <- X_dot_i[ind_RE_i, ind_FE_i]
    b_i <- b_mat[i, ]
    u_i <- b_i[ind_RE_i] + X_dot_i %*% betas_vec[ind_FE_i];
    D_inv_i <- D_inv[[patt_i]]
    XD_i <- t(X_dot_i) %*% D_inv_i
    XDX_i <- XD_i %*% X_dot_i
    J_i <- diag(p_HC)[, ind_FE_i]
    sum_JXDu  = sum_JXDu  + J_i %*% (XD_i %*% u_i)
    sum_JXDXJ = sum_JXDXJ +  J_i %*% XDX_i %*% t(J_i)
  }

  Sigma_1  = solve(prior_Tau_betas_HC + sum_JXDXJ)
  mean_1 = Sigma_1 %*% (prior_Tau_betas_HC %*% prior_mean_betas_HC + sum_JXDu) 
  
  list(D_inv = D_inv, mean_1 = mean_1, Sigma_1 = Sigma_1)
  
}

sds <- sqrt(diag(D))
R <- cov2cor(D)
L <- chol(R)

U <- chol(D)
b_mat <- do.call(cbind, b)
prior_mean_betas_HC <- mean_betas_HC
prior_Tau_betas_HC <- Tau_betas_HC

testR <- update_betasR(betas, 
                      prior_mean_betas_HC, prior_Tau_betas_HC,
                      b_mat, D, X_dot,
                      ind_FE,
                      ind_FE_HC,
                      id_patt,
                      ind_RE_patt,
                      ind_FE_patt,
                      it = 1,
                      y)

testR$mean_1
testR$Sigma_1
testR$D_inv
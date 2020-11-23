create_X_dot3 <- function(nres, nfes_HC, z_in_x, x_in_z, X_HC, nT, unq_idL, xbas_in_z) {
  
  n_outcomes <- length(nres) # number of outcomes
  n_res <- sum(nres) # total number of RE
  
  M <- matrix(0, nrow= n_res*nT, ncol= sum(nfes_HC))
  
  for (j in seq_len(n_outcomes)) { # j-th outcome
    
    ids <- unq_idL[[j]] # ids present in outcome-j
    ids_rows <- (ids-1) * n_res # 1st row for each id
    
    rows <- sum(nres[1:j-1]) + z_in_x[[j]] + rep(ids_rows, each= length(z_in_x[[j]]))
    cols <- rep(sum(nfes_HC[1:j-1]) + x_in_z[[j]], times= length(ids))
    M[cbind(rows, cols)] <- 1 # add 1 for each z_in_x
    
    bas_cols <- xbas_in_z[[j]]
    
    for (k in z_in_x[[j]]) { # k-th RE in z_in_x
      
      if(sum(bas_cols[k, ])==0) next
      
      M[sum(nres[1:j-1]) + k + ids_rows, sum(nfes_HC[1:j-1]) + which(bas_cols[k, ])] <- X_HC[[j]][[k]]
    }
  }
  M
}


X_dot <- create_X_dot3(nres, nfes_HC, z_in_x, x_in_z, X_HC, nT, unq_idL, xbas_in_z)

#

set.seed(2021)
id <- sample(seq_len(nT), 1)
id <- 2
rows <- lapply(idL, match, x=id) 
mapply(function(XX, r){XX[r,]}, X, rows)

X_dot[seq_len(sum(nres)) + (id-1)*sum(nres),]

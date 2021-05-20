JMbayes2:::linpred_surv(X_H, betas, Z_H, b, id_H)
X <- X_H
betas <- betas
Z <- Z_H
b <- b
id <- id_H
i <- 1
function (X, betas, Z, b, id) {
  out <- vector("list", length(X))
  for (i in seq_along(X)) {
    X_i <- X[[i]]
    Z_i <- Z[[i]]
    betas_i <- betas[[i]]
    b_i <- b[[i]]
    id_i <- id[[i]]
    out[[i]] <- matrix(0.0, nrow = nrow(X_i[[1]]), ncol = length(X_i))
    for (j in seq_along(X_i)) {
      X_ij <- X_i[[j]]
      Z_ij <- Z_i[[j]]
      out[[i]][, j] <- X_ij %*% betas_i + rowSums(Z_ij * b_i[id_i, ])
    }
  }
  out
}


unique(id_i)
b_i[1, ]

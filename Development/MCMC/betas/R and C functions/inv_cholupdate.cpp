#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat rank1_update (const mat &U, // performs rank-1 update. If U = chol(M), returns chol(M + v * v.t())
                  const vec &v) {
  
  uword n = v.n_elem;
  mat Res = U;
  vec v2  = v;
  
  for (uword i = 0; i < n; i++) {
    double r = pow( pow(Res.at(i, i), 2) + pow(v2.at(i), 2), 0.5);
    double c = r / Res.at(i, i);
    double s = v2.at(i) / Res.at(i, i);
    Res.at(i, i) = r;
    
    if (i < n-1) {
      Res.submat(i, i + 1, i, n - 1) = (Res.submat(i, i + 1, i, n - 1) + s * v2.rows(i + 1, n - 1).t()) / c;
      v2.rows(i + 1, n - 1) = c * v2.rows(i + 1, n - 1) - s * Res.submat(i, i + 1, i, n - 1).t();
    }
  }
  return Res;
}

// [[Rcpp::export]]
mat rank1_updateL(const mat &L, // performs rank-1 update. If L = chol(M).t(), returns chol(M + v * v.t()).t()
                 const vec &v) {
  
  uword n = v.n_elem;
  mat Res = L;
  vec v2  = v;
  
  for (uword i = 0; i < n; i++) {
    double r = pow( pow(Res.at(i, i), 2) + pow(v2.at(i), 2), 0.5);
    double c = r / Res.at(i, i);
    double s = v2.at(i) / Res.at(i, i);
    Res.at(i, i) = r;
    
    if (i < n-1) {
      Res.submat(i + 1, i, n - 1, i) = (Res.submat(i + 1, i, n - 1, i) + s * v2.rows(i + 1, n - 1)) / c;
      v2.rows(i + 1, n - 1) = c * v2.rows(i + 1, n - 1) - s * Res.submat(i + 1, i, n - 1, i);
    }
  }
  return Res;
}

// [[Rcpp::export]]
mat chol_update(const mat &U, // If U = chol(M), returns chol(M.submat(keep, keep))
                const uvec &keep) { // keep must be a sorted vector, i.e, {2, 4, 5}, and counts from 0
  
  // later we can try to extend this approach further to obtain inv(U) from the required inv(U_i)
  
  uvec rem = regspace<uvec>(0,  U.n_cols - 1); rem.shed_rows(keep); // cols-rows to remove
  mat Res = U;
  uword n = rem.n_elem;
  
  for (uword i = 0; i < n; i++) { // rank-1 update for each col-row to be removed
    
    uword last_col = Res.n_cols - 1;
    
    if(rem.at(i) < last_col) {
      Res.submat(rem.at(i) + 1, rem.at(i) + 1, last_col, last_col) = rank1_update(Res.submat(rem.at(i) + 1, rem.at(i) + 1, last_col, last_col),
                 Res.submat(rem.at(i), rem.at(i) + 1, rem.at(i), last_col).t());
    }
    
    Res.shed_row(rem.at(i));
    Res.shed_col(rem.at(i));
    rem = rem - 1;
  }
  return Res;
} 

// [[Rcpp::export]]
mat inv_trimatu (mat U) {
  
  return inv(trimatu(U));
  
} 


/*** R
set.seed(2020)
nres <- 6
R <- matrix(0, sum(nres), sum(nres))
R[lower.tri(R, diag= FALSE)] <- runif( (sum(nres)^2 - sum(nres))/2, 0.2, 0.8)
R <- R + t(R)
diag(R) <- 1
sd <- runif(nres, 0.5, 3)
D <- diag(sd) %*% R %*% diag(sd)
U <- chol(D)
L <- t(U)

# rank-1 update
v <- runif(nres) # random vector
t(rank1_updateL(L, v))
t(ramcmc::chol_update(L, v))
rank1_update(U, v)
chol(D + v %*% t(v))

# chol_update 1 
keep <- c(1, 3, 5)
chol_update(U, keep-1)
chol(D[keep, keep])

# inv_trimatu
inv_trimatu(chol_update(U, keep-1))
solve(chol(D[keep, keep]))

# chol_update, all possible ind combinations
all_keep <- lapply(seq_len(nres), combn, x = seq_len(nres), simplify = FALSE)
all_keep <- unlist(all_keep, recursive = FALSE)

for(keep in all_keep) {
  
  C <- chol_update(U, keep-1)
  R <- chol(D[keep, keep])
  
  cat(paste0(all.equal(C, R), " "))
  if(!all.equal(C, R)) print(paste0(keep, "\n"))
  
}

*/

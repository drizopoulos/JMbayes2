#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

mat rnorm_mat (const uword& rows, const uword& cols) {
  mat out(rows, cols);
  out.each_col([&](vec& x) {x = as<vec>(rnorm(rows)); } );
  return out;
}

// [[Rcpp::export]]
mat propose_mvnorm_vec (const int &n, // number of samples
                        const mat &U, // upper triangular Choleski factorization of the vcov matrix
                        const double &scale) {
  uword ncols = U.n_cols;
  mat res(n, ncols);
  res = scale * (rnorm_mat(n, ncols) * U);
  return res.t(); // j-th column reports the j-th sample
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

sample <- propose_mvnorm_vec(1000000, U, 1)
cov(t(sample))
D
solve(D)

rowMeans(sample)
*/

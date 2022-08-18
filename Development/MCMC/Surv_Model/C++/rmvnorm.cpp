#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
vec propose_mvnorm_vec (const mat &Sigma) {
    mat U = chol(Sigma, "lower");
    vec res = U * randn(U.n_cols);
    return res;
}

/*** R
SS <- cbind(c(5, 1.5), c(1.5, 3))
M <- 150000
mat <- matrix(0.0, M, 2)
for (m in seq_len(M)) {
    mat[m, ] <- c(propose_mvnorm_vec(SS))
}
var(mat)
*/


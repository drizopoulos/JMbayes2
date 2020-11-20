//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
cube foo (const mat& x) {
  cube out(x.n_cols, x.n_cols, x.n_rows);
  for (uword i = 0; i < x.n_rows; i++) {
    out.slice(i) = kron(x.row(i), x.row(i).t());
  }
  return out;
}

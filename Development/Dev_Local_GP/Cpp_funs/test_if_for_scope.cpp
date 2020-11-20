//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

field<vec> foo (field<vec> &X, bool t) {
  if (t) {
    uword n = X.n_elem;
    uword k = X.at(0).n_elem;
    for (uword i = 0; i < k; i++) {
      for (uword j = 0; j < n; j++) {
        X.at(j).at(i) = X.at(j).at(i) - 10;
      }
    }
  }
  return(X);
}
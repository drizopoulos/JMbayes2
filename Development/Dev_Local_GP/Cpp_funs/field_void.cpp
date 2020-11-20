//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void foo (field<arma::vec>& X) {
  uword n = X.n_elem;
  for (uword i = 0; i < n; i++) {
    X.at(i) = X.at(i) + 1;
  }
}

// [[Rcpp::export]]
field<arma::vec> test_void (field<arma::vec>& X) {
  foo(X);
  return X;
}
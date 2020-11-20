//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


void defoo (field<arma::vec>& Y) {
  uword n = Y.n_elem;
  for (uword i = 0; i < n; i++) {
    Y.at(i) = Y.at(i) + 10;
  }
}

void foo (field<arma::vec>& X) {
  uword n = X.n_elem;
  defoo(X);
  for (uword i = 0; i < n; i++) {
    X.at(i) = X.at(i) + 1;
  }
}

// [[Rcpp::export]]
field<arma::vec> test_void (field<arma::vec>& X) {
  foo(X);
  return X;
}
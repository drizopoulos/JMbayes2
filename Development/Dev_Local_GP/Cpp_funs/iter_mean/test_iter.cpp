//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat foo (field<mat> x) {
  mat out(size(x.at(0)));
  uword n = x.n_elem;
  for (uword i = 1; i < n; i++) {
    out = x.at(i) - out;
  }
  return(out);
}
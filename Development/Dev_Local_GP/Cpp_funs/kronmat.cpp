//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat kronmat (const vec& x) {
  mat out = kron(x, x.t());
  return out;
}

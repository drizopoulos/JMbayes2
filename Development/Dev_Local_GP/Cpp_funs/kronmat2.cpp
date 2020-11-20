//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat kronmat2 (const vec& x) {
  mat out = x * x.t();
  return out;
}

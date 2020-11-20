//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

cube foo2 () {
  cube x(500, 4, 1, fill::zeros);
  x.set_size(500, 4, 100001);
  return x;
}
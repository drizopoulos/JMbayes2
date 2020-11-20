//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]

cube foo () {
  cube x(500, 4, 1, fill::zeros);
  x.insert_slices(0, 100000);
  return x;
}
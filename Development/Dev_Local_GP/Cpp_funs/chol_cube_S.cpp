//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
cube chol_cube (const cube &S) {
  cube out = S;
  out.each_slice([](mat &X){X = chol(X, "lower");});
  return out;
}
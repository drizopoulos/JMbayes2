//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
cube choleach (const cube& S) {
  cube S_chol = S;
  S_chol.each_slice( [&] (mat& X) {chol(X); } );
  return S_chol;
}
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
cube cholfor (const int& n, const cube& S) {
  uword ncol_per_slice = S.n_cols;
  uword slices = S.n_slices;
  cube out(S);
  for (uword i = 0; i < slices; i++) {
    out.slice(i) =  chol(S.slice(i));
  }
  return out;
}
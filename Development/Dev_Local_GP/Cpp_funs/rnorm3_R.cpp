//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
mat rnorm3 (const int& rows, const int& cols) {
  mat out(rows, cols, fill::randn);
  return out;
}

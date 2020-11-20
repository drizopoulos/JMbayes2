//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
cube res (bool x) {
  cube out();
  if (x) {
    out(493, 4, 3000, fill::zeros);
    return out;
  } else {
    cube out(493, 4, 1, fill::zeros);
    return(out);
  }
  return out;
}
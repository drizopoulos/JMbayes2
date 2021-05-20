//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
int check (int x, int y) {
  if (x > y - 1) {
    for (int i = 0; i < 10; i++) {
      x += i;
    }
  }
  return x;
}
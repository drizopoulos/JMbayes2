//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List foo (mat &res) {
  cube res_b(10, 4, 20, fill::zeros);
  //res_b.slice(0) = res;
  //return res_b.slices(0, 2 - 1);
  cube out = res_b.slices(5, 9);
  return List::create(
    Named("b") = out
  );
}
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
mat rnorm2 (const int& rows, const int& cols, const double& mu, const double& sigma) {
  mat out(rows, cols);
  out.each_col([&](vec& x) {x = as<vec>(rnorm(rows)); } );
    return out;
}
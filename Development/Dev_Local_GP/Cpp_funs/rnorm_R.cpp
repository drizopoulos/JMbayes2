//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
mat rnorm (const int& rows, const int& cols, const double& mu, const double& sigma) {
  mat out(rows, cols);
  out.imbue([&]() {return R::rnorm(mu, sigma); } );
  //out.each_col([&](vec& x) {x = rnorm(n); } );
  return out;
}

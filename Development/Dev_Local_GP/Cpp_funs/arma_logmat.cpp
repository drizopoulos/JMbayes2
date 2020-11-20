//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double arma_logmat (const mat& M, const uword &row, const uword &col) {
  return std::log(M.at(row, col));
}

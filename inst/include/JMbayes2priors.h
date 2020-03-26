#ifndef JMBAYES2PRIORS_H
#define JMBAYES2PRIORS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

arma::vec dlkj(arma::mat Omega, arma::vec eta = 1, bool logscale = true) {
  arma::vec out;
  if (logscale = true) {
    out = pow(det(Omega), eta - 1);
  } else {
    out = exp(pow(det(Omega), eta - 1);)
  }
  return(out);
}

arma::mat rwishart(int df, const arma::mat& S) {
  int m = S.n_rows;
  mat Z(m, m);
  for (int i = 0; i < m; ++i) {
    Z(i, i) = sqrt(R::rchisq(df - i));
  }
  for (int j = 0; j < m; ++j){
    for (int i = j + 1; i < m; ++i){
      Z(i, j) = R::rnorm(0, 1);
    }
  }
  arma::mat C = arma::trimatl(Z).t() * arma::chol(S);
  return C.t() * C;
}

#endif
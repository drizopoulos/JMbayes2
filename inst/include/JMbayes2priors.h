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

#endif
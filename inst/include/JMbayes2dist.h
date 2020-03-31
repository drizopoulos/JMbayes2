#ifndef JMBAYES2PRIORS_H
#define JMBAYES2PRIORS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// LKJ
arma::vec dlkj(arma::mat Omega, arma::vec eta = 1, bool logscale = true) {
  arma::vec out;
  if (logscale = true) {
    out = pow(det(Omega), eta - 1);
  } else {
    out = exp(pow(det(Omega), eta - 1);)
  }
  return(out);
}

// symmetric Dirichlet
arma::vec ddirichlet_sym(const arma::vec& x, const arma::vec alpha, bool logscale = true) {
  //int K = alpha.n_elem;
  arma::vec out;
  if (logscale == true) {
    out = lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) % log(x));
  } else {
    out = exp(lgamma(sum(alpha)) - sum(lgamma(alpha)) + sum((alpha - 1) % log(x)));
  }
  return(out);
}

#endif

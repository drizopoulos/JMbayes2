#include "../inst/include/JMbayes2cpp.h"
#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List mcmc_jm2 (List initials, List Data, List priors, List scales, List Covs,
               List control, bool interval_cens, bool multiState) {
  
}
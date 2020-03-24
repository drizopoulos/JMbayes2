#include "../inst/include/JMbayes2cpp.h"
#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
int mcmc_jm2 (List initials, List Data, List priors, List scales, List Covs,
               List control, bool interval_cens, bool multiState) {
  
  // Extract Data
  List idL = as<List>(Data["idL"]);
  field<uvec> idLF = List2Field_vec(idL);
  
  
  
  return(0)
}
#include <RcppArmadillo.h>
#include "../inst/include/JMbayes2cpp.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
int mcmc_jm2 (List initials, List Data, List priors, List scales, List Covs,
               List control, bool interval_cens, bool multiState) {
  
  // Extract Data
  List idL = as<List>(Data["idL"]);
  field<uvec> idL_F = List2Field_uvec(idL);
  List idL_lp = as<List>(Data["idL_lp"]);
  field<uvec> idL_lp_F = List2Field_uvec(idL_lp);
  List unq_idL = as<List>(Data["unq_idL"]);
  field<uvec> unq_idL_F = List2Field_uvec(unq_IdL);
  List y = as<List>(Data["y"]);
  field<vec> y_F = List2Field_vec(y);
  List X = as<List>(Data["X"]);
  field<mat> X_F = List2Field_mat(X);
  List Z = as<List>(Data["Z"]);
  field<mat> Z_F = List2Field_mat(Z);
  List Xhc = as<List>(Data["Xhc"]);
  field<mat> Xhc_F = List2Field_mat(Xhc);
  List columns_HC = as<List>(Data["columns_HC"]);
  field<uvec> columns_HC_F = List2Field_uvec(columns_HC);
  // columns_nHC can be NULL
  //List columns_nHC = as<List>(Data["columns_nHC"]);
  //field<uvec> columns_nHC_F = List2Field_uvec(columns_nHC);
  // idT currently factor/character from R
  List Time_right = as<List>(Data["Time_right"]);
  field<vec> Time_right_F = List2Field_vec(Time_right);
  List Time_left = as<List>(Data["Time_left"]);
  field<vec> Time_left_F = List2Field_vec(Time_left);
  
  
  
  
  return(0)
}
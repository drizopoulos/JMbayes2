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
  // CharacterVector idT = as<CharacterVector>(Data["idT"]);
  List Time_right = as<List>(Data["Time_right"]);
  field<vec> Time_right_F = List2Field_vec(Time_right);
  List Time_left = as<List>(Data["Time_left"]);
  field<vec> Time_left_F = List2Field_vec(Time_left);
  uvec delta = as<uvec>(Data["delta"]);
  uvec which_event = as<uvec>(Data["which_event"]);
  uvec which_right = as<uvec>(Data["which_right"]);
  uvec which_left = as<uvec>(Data["which_left"]);
  uvec which_interval = as<uvec>(Data["which_interval"]);
  mat W0_H = as<mat>(Data["W0_H"]);
  mat W_H = as<mat>(Data["W_H"]);
  List X_H = as<List>(Data["X_H"]);
  field<mat> X_H_F = List2Field_mat(X_H);
  List Z_H = as<List>(Data["Z_H"]);
  field<mat> Z_H_F = List2Field_mat(Z_H);
  List U_H = as<List>(Data["U_H"]);
  field<mat> U_H_F = List2Field_mat(U_H);
  mat W0_h = as<mat>(Data["W0_h"]);
  mat W_h = as<mat>(Data["W_h"]);
  List X_h = as<List>(Data["X_h"]);
  field<mat> X_h_F = List2Field_mat(X_h);
  List Z_h = as<List>(Data["Z_h"]);
  field<mat> Z_h_F = List2Field_mat(Z_h);
  List U_h = as<List>(Data["U_h"]);
  field<mat> U_h_F = List2Field_mat(U_h);
  mat W0_H2 = as<mat>(Data["W0_H2"]);
  mat W_H2 = as<mat>(Data["W_H2"]);
  List X_H2 = as<List>(Data["X_H2"]);
  field<mat> X_H2_F = List2Field_mat(X_H2);
  List Z_H2 = as<List>(Data["Z_H2"]);
  field<mat> Z_H2_F = List2Field_mat(Z_H2);
  List U_H2 = as<List>(Data["U_H2"]);
  field<mat> U_H2_F = List2Field_mat(U_H2);
  vec log_Pwk = as<vec>(Data["log_Pwk"]);
  vec log_Pwk2 = as<vec>(Data["log_Pwk"]);
  
  
  return(0)
}
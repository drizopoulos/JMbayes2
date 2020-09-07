#ifndef JMBAYES2SURV_H
#define JMBAYES2SURV_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_b (cube &b, mat &b_mat, const field<mat> &Xbetas, const field<mat> &Z, const field<uvec> &id, 
               const field<mat> &y, const vec &scales, const vec &extra_parms, 
               const CharacterVector &families, const CharacterVector &links, const field<uvec> &ids,
               const field<uvec> &unq_ids, const mat &L, , 
               const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
               const vec &W0H2_bs_gammas, const vec &WH_gammas,
               const vec &Wh_gammas, const vec &WH2_gammas,
               const vec &WlongH_alphas, const vec &Wlongh_alphas,
               const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
               const uvec &indFast_H, const uvec &which_event,
               const uvec &which_right_event, const uvec &which_left,
               const bool &any_interval, const uvec &which_interval, 
               const field<uvec> ind_RE) {
  mat proposed_b = propose_mvnorm_mat(1, L, scales) + b_mat;
  field<mat> proposed_b_field = mat2field_mat(proposed_b, ind_RE); 
  vec numerator_b = log_b(Xbetas, Z, proposed_b_field, proposed_b,
                          id, y, scales,
                          extra_parms, families,
                          links, ids,
                          unq_ids, L, 
                          W0H_bs_gammas, W0h_bs_gammas,
                          W0H2_bs_gammas, WH_gammas,
                          Wh_gammas, WH2_gammas,
                          WlongH_alphas, Wlongh_alphas,
                          WlongH2_alphas, log_Pwk, log_Pwk2,
                          indFast_H, which_event,
                          which_right_event, which_left,
                          any_interval, which_interval); 
}
#ifndef JMBAYES2SURV_H
#define JMBAYES2SURV_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

vec log_b (const field<mat> &X, const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2, 
           const field<mat> &Z, const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2, 
           const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2, 
           const field<vec> &betas,
           const field<mat> &b, const mat &b_mat, 
           const field<mat> &y, const field<vec> &eta, const vec &sigmas,
           const vec &extra_parms, const CharacterVector &families,
           const CharacterVector &links, const field<uvec> &ids,
           const field<uvec> &unq_ids,
           const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
           const vec &W0H2_bs_gammas, const vec &WH_gammas,
           const vec &Wh_gammas, const vec &WH2_gammas,
           const vec &WlongH_alphas, const vec &Wlongh_alphas,
           const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
           const uvec &indFast_H, const uvec &which_event,
           const uvec &which_right_event, const uvec &which_left,
           const bool &any_interval, const uvec &which_interval,
           const mat &L, const vec &sds) {
  // log_lik_y
  field<vec> eta = linpred_mixed(X, betas, Z, b, id);
  vec logLik_long = log_long(y, eta, sigmas, extra_parms, families, links, ids, unq_ids);
  // logLik_re
  vec logLik_re = log_re(b_mat, L, sds);
  // loglik_surv
  field<mat> eta_H = linpred_surv(X_H, betas, Z_H, b, id_H);
  Wlong_H = docall_cbindF(create_Wlong(eta_H, FunForms, U_H, FunForms_ind));
  if (any_event) {
    field<mat> eta_h = linpred_surv(X_h, betas, Z_h, b, id_h);
    Wlong_h = docall_cbindF(create_Wlong(eta_h, FunForms, U_h, FunForms_ind));
  }
  if (any_interval) {
    field<mat> eta_H2 = linpred_surv(X_H2, betas, Z_H2, b, id_H);
    Wlong_H2 = docall_cbindF(create_Wlong(eta_H2, FunForms, U_H2, FunForms_ind));
  }
  vec WlongH_alphas = Wlong_H * alphas;
  vec Wlongh_alphas(W0_h.n_rows);
  vec WlongH2_alphas(W0_H2.n_rows);
  if (any_event) {
    Wlongh_alphas = Wlong_h * alphas;
  }
  if (any_interval) {
    WlongH2_alphas = Wlong_H2 * alphas;
  }
  vec log_lik_surv = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                              WH_gammas, Wh_gammas, WH2_gammas,
                              WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                              log_Pwk, log_Pwk2, indFast_H,
                              which_event, which_right_event, which_left,
                              any_interval, which_interval);
  vec out = logLik_long + log_re + log_lik_surv;
  return out;
}

/*
 vec log_b (const field<mat> &X, const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
 const field<mat> &Z, const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
 const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
 const field<vec> &betas,
 const field<mat> &b, const mat &b_mat,
 const field<mat> &y, const vec &scales, const mat &chol_S,
 const vec &extra_parms, const CharacterVector &families, const CharacterVector &links,
 const field<uvec> &id, const uvec &id_H, const uvec &id_h,
 const field<uvec> &ids, const field<uvec> &unq_ids,
 const mat &W0_h, const mat &W0_H2,
 const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
 const vec &W0H2_bs_gammas, const vec &WH_gammas,
 const vec &Wh_gammas, const vec &WH2_gammas,
 const vec &alphas,
 const vec &log_Pwk, const vec &log_Pwk2,
 const uvec &indFast_H, const uvec &which_event,
 const uvec &which_right_event, const uvec &which_left,
 const bool &any_interval, const uvec &which_interval,
 const field<uvec> &FunForms, const field<uvec> &FunForms_ind) {
 // log_lik_y
 field<vec> eta = linpred_mixed(X, betas, Z, b, id);
 vec logLik_long = log_long(y, eta, scales, extra_parms, families, links, ids, unq_ids);
 // log_pb
 vec log_re = log_dmvnrm_chol(b_mat, chol_S);
 // log_lik_surv
 field<mat> eta_H = linpred_surv(X_H, betas, Z_H, b, id_H);
 Wlong_H = docall_cbindF(create_Wlong(eta_H, FunForms, U_H, FunForms_ind));
 if (any_event) {
 field<mat> eta_h = linpred_surv(X_h, betas, Z_h, b, id_h);
 Wlong_h = docall_cbindF(create_Wlong(eta_h, FunForms, U_h, FunForms_ind));
 }
 if (any_interval) {
 field<mat> eta_H2 = linpred_surv(X_H2, betas, Z_H2, b, id_H);
 Wlong_H2 = docall_cbindF(create_Wlong(eta_H2, FunForms, U_H2, FunForms_ind));
 }
 vec WlongH_alphas = Wlong_H * alphas;
 vec Wlongh_alphas(W0_h.n_rows);
 vec WlongH2_alphas(W0_H2.n_rows);
 if (any_event) {
 Wlongh_alphas = Wlong_h * alphas;
 }
 if (any_interval) {
 WlongH2_alphas = Wlong_H2 * alphas;
 }
 vec log_lik_surv = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
 WH_gammas, Wh_gammas, WH2_gammas,
 WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
 log_Pwk, log_Pwk2, indFast_H,
 which_event, which_right_event, which_left,
 any_interval, which_interval);
 //
 vec out = logLik_long + log_re + log_lik_surv;
 return out;
 }
 */



void update_b (field<mat> &b, mat &b_mat, const cube &chol_S, const vec &scale_sigmas, ind_RE) {
  
  mat proposed_b = propose_mvnorm_mat(1, L, scale_b) + b_mat;
  field<mat> proposed_b_field = mat2field_mat(proposed_b, ind_RE); 
  
  
  
  vec numerator_b = log_b(Xbetas, Z, proposed_b_field, proposed_b,
                          id, y, scale_b,
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
  vec denominator_b = log_b(Xbetas, Z, b, b_mat,
                            id, y, scale_b,
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
  vec log_ratio = numerator_b - denominator_b;
  uword n = log_ratio.n_elem;
  for (uword i = 0; i < n; i++) {
    if (std::isfinite(log_ratio.at(i)) && exp(log_ratio.at(i)) > R::runif(0, 1)) {
      acceptance_b.at(i, it) = 1;
      b_mat.row(i) = proposed_b.row(i);
    }
    if (it > 19) {
      scale_b.at(i) =
        robbins_monro(scale_b.at(i),
                      acceptance_b.at(i, it), it);
    }
    res_b.row(i) = b_mat.row(i);
  }
  b = mat2field_mat(b_mat, ind_RE);
}
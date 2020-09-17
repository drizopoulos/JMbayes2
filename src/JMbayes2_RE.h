#ifndef JMBAYES2SURV_H
#define JMBAYES2SURV_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

vec log_b (const field<mat> &X, const field<mat> &Z, const field<vec> &betas, 
           const field<mat> &b, const mat &b_mat, 
           const field<uvec> &id, const field<mat> &y, const vec &sigmas,
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
  // log_lik_long
  field<vec> eta = linpred_mixed(X, betas, Z, b, id);
  vec logLik_long = log_long(y, eta, sigmas, extra_parms, families, links, ids, unq_ids);
  // loglik_surv
  vec logLik_surv = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                             WH_gammas, Wh_gammas, WH2_gammas,
                             WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                             log_Pwk, log_Pwk2, indFast_H,
                             which_event, which_right_event, which_left,
                             any_interval, which_interval);
  // logLik_re
  vec logLik_re = log_re(b_mat, L, sds);
  // out
  vec out = logLik_long + logLik_re + logLik_surv;
  return out;
}

void update_b (field<mat> &b, mat &b_mat, const cube &chol_S, const vec &scale_b, 
               const field<uvec> &ind_RE, mat &Wlong_H, mat &Wlong_h, mat &Wlong_H2,
               const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
               const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
               const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
               const field<vec> &betas, const vec &alphas, 
               const uvec &id_H, const uvec &id_h,
               const field<uvec> &FunForms, const field<uvec> &FunForms_ind,
               const bool &any_event, const bool &any_interval, 
               const field<mat> &X, const field<mat> &Z, 
               const field<uvec> &id, const field<mat> &y,  const vec &sigmas,
               const vec &extra_parms, const CharacterVector &families,
               const CharacterVector &links, const field<uvec> &ids,
               const field<uvec> &unq_ids, const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
               const vec &W0H2_bs_gammas, const vec &WH_gammas,
               const vec &Wh_gammas, const vec &WH2_gammas,
               const vec &log_Pwk, const vec &log_Pwk2,
               const uvec &indFast_H, const uvec &which_event,
               const uvec &which_right_event, const uvec &which_left,
               const bool &any_interval, const uvec &which_interval,
               const mat &L, const vec &sds, vec &denominator_b, 
               const uword &it, const uword &n_rows_W0_h, 
               const uword &n_rows_W0_H2, 
               ) {
  // propose new random effects
  mat proposed_b = propose_mvnorm_mat(1, chol_S, scale_b) + b_mat;
  field<mat> proposed_b_field = mat2field_mat(proposed_b, ind_RE); 
  update_Wlong(Wlong_H, Wlong_h, Wlong_H2, X_H, X_h, X_H2, Z_H, Z_h, Z_H2, 
               U_H, U_h, U_H2, betas, proposed_b_field, id_H, id_h, FunForms, FunForms_ind, 
               any_event, any_interval);
  vec WlongH_alphas_updated = Wlong_H * alphas;
  vec Wlongh_alphas_updated(n_rows_W0_h);
  vec WlongH2_alphas_updated(n_rows_W0_H2);
  if (any_event) {
    Wlongh_alphas_updated = Wlong_h * alphas;
  }
  if (any_interval) {
    WlongH2_alphas_updated = Wlong_H2 * alphas;
  }
  vec numerator_b = log_b(X, Z, betas, proposed_b_field, proposed_b, id, y, sigmas, 
                          extra_parms, families, links, ids, unq_ids, W0H_bs_gammas, 
                          W0h_bs_gammas, W0H2_bs_gammas, WH_gammas, 
                          Wh_gammas, WH2_gammas, WlongH_alphas_updated, 
                          Wlongh_alphas_updated, WlongH2_alphas_updated, 
                          log_Pwk, log_Pwk2, indFast_H, which_event, which_right_event, 
                          which_left, any_interval, which_interval, 
                          L, sds);
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
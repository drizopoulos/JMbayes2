#ifndef JMBAYES2RE_H
#define JMBAYES2RE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
#include "JMbayes2_Long.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_b (field<mat> &b, mat &b_mat, field<vec> &eta,
               vec &logLik_long, vec &logLik_surv, vec &logLik_re, 
               mat &Wlong_H, mat &Wlong_h, mat &Wlong_H2, 
               const cube &chol_S, vec &scale_b, 
               const field<uvec> &ind_RE, 
               const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
               const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
               const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
               const field<vec> &betas, const vec &alphas, 
               const uvec &id_H, const uvec &id_h,
               const field<uvec> &FunForms, const field<uvec> &FunForms_ind,
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
               const uvec &which_interval, const bool &any_event, const bool &any_interval, 
               const mat &L, const vec &sds,
               const uword &it, const uword &n_rows_W0_h, 
               const uword &n_rows_W0_H2, const field<uvec> &rows_Wlong_H,
               mat &acceptance_b, cube &res_b
               ) {
  // calculate denominator_b
  vec denominator_b = logLik_long + logLik_surv + logLik_re;
  // propose new random effects in mat and field<mat> form
  mat proposed_b = propose_mvnorm_mat(1, chol_S, scale_b) + b_mat;
  field<mat> proposed_b_field = mat2field_mat(proposed_b, ind_RE);
  // calculate log_lik_long based on proposed_b
  field<vec> eta_proposed = linpred_mixed(X, betas, Z, proposed_b_field, id);
  vec logLik_long_proposed = log_long(y, eta_proposed, sigmas, extra_parms, families, links, ids, unq_ids);
  // create copies of Wlong_H, Wlong_h and Wlong_H2 to be updated
  mat Wlong_H_proposed(size(Wlong_H));
  mat Wlong_h_proposed(size(Wlong_h));
  mat Wlong_H2_proposed(size(Wlong_H2));
  // update Wlong_H_updated, Wlong_h_updated and Wlong_H2_updated based on proposed_b
  update_Wlong(Wlong_H_proposed, Wlong_h_proposed, Wlong_H2_proposed, X_H, X_h, X_H2, Z_H, Z_h, Z_H2, 
               U_H, U_h, U_H2, betas, proposed_b_field, id_H, id_h, FunForms, FunForms_ind, 
               any_event, any_interval);
  vec WlongH_alphas_proposed = Wlong_H_proposed * alphas;
  vec Wlongh_alphas_proposed(n_rows_W0_h);
  vec WlongH2_alphas_proposed(n_rows_W0_H2);
  if (any_event) {
    Wlongh_alphas_proposed = Wlong_h_proposed * alphas;
  }
  if (any_interval) {
    WlongH2_alphas_proposed = Wlong_H2_proposed * alphas;
  }
  // calculate logLik_Surv_proposed
  vec logLik_surv_proposed = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas, 
                                      WH_gammas, Wh_gammas, WH2_gammas, 
                                      WlongH_alphas_proposed, Wlongh_alphas_proposed, 
                                      WlongH2_alphas_proposed, 
                                      log_Pwk, log_Pwk2, indFast_H, 
                                      which_event, which_right_event, which_left, 
                                      any_interval, which_interval);
  // logLik_re
  vec logLik_re_proposed = log_re(b_mat, L, sds);
  // calculate the numerator
  vec numerator_b = logLik_long_proposed + logLik_surv_proposed + logLik_re_proposed;
  // log_rati
  vec log_ratio = numerator_b - denominator_b;
  uword n = log_ratio.n_elem;
  for (uword i = 0; i < n; i++) {
    if (std::isfinite(log_ratio.at(i)) && exp(log_ratio.at(i)) > R::runif(0, 1)) {
      acceptance_b.at(it, i) = 1;
      b_mat.row(i) = proposed_b.row(i);
      logLik_long.at(i) = logLik_long_proposed.at(i);
      logLik_surv.at(i) = logLik_surv_proposed.at(i);
      logLik_re.at(i) = logLik_re_proposed.at(i);
      uword n_outcomes = eta.n_elem;
      for (uword j = 0; j < n_outcomes; j++) {
        eta.at(j).at(i) = eta_proposed.at(j).at(i);
      }
      Wlong_H.rows(rows_Wlong_H.at(i)) = Wlong_H_proposed.rows(rows_Wlong_H.at(i));
      Wlong_h.row(i) = Wlong_h_proposed.row(i);
      Wlong_H2.row(i) = Wlong_H2_proposed.row(i);
    }
    if (it > 19) {
      scale_b.at(i) =
        robbins_monro(scale_b.at(i),
                      acceptance_b.at(i, it), it);
    }
  }
  res_b.slice(it) = b_mat;
  b = mat2field_mat(b_mat, ind_RE);
}

#endif
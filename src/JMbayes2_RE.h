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
               vec &WlongH_alphas, vec &Wlongh_alphas, vec &WlongH2_alphas,
               const cube &chol_S, vec &scale_b,
               const field<uvec> &ind_RE,
               const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
               const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
               const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
               const mat &Wlong_bar, const field<vec> &betas, const vec &alphas,
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
               const uvec &which_interval, const bool &any_event,
               const bool &any_interval,
               const mat &L, const vec &sds,
               const uword &it, const field<uvec> &rows_Wlong_H,
               const field<uvec> &idL_ind,
               mat &acceptance_b, cube &res_b
               ) {
  // calculate denominator_b
  vec denominator_b = logLik_long + logLik_surv + logLik_re;
  // propose new random effects in mat and field<mat> form
  mat proposed_b_mat = propose_mvnorm_mat(1, chol_S, scale_b) + b_mat;
  field<mat> proposed_b = mat2field_mat(proposed_b_mat, ind_RE);
  // calculate log_lik_long based on proposed_b_mat
  field<vec> eta_proposed = linpred_mixed(X, betas, Z, proposed_b, id);
  vec logLik_long_proposed = log_long(y, eta_proposed, sigmas, extra_parms,
                                      families, links, ids, unq_ids);
  // calculate Wlong_H, Wlong_h and Wlong_H2 based on the proposed_b
  mat Wlong_H_proposed =
    calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, betas, proposed_b, id_H, FunForms,
                    FunForms_ind);
  mat Wlong_h_proposed =
    calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, betas, proposed_b, id_h, FunForms,
                    FunForms_ind);
  mat Wlong_H2_proposed =
    calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, betas, proposed_b, id_H,
                    FunForms, FunForms_ind);
  // create and initiate  WlongH_alphas_proposed, Wlongh_alphas_proposed, WlongH2_alphas_proposed
  vec WlongH_alphas_proposed = Wlong_H_proposed * alphas;
  vec Wlongh_alphas_proposed(Wlongh_alphas.n_rows);
  vec WlongH2_alphas_proposed(WlongH2_alphas.n_rows);
  if (any_event) {
    Wlongh_alphas_proposed = Wlong_h_proposed * alphas;
  }
  if (any_interval) {
    WlongH2_alphas_proposed = Wlong_H2_proposed * alphas;
  }
  // calculate logLik_Surv_proposed
  vec logLik_surv_proposed =
    log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
             WH_gammas, Wh_gammas, WH2_gammas,
             WlongH_alphas_proposed, Wlongh_alphas_proposed,
             WlongH2_alphas_proposed,
             log_Pwk, log_Pwk2, indFast_H,
             which_event, which_right_event, which_left,
             any_interval, which_interval);
  // logLik_re
  vec logLik_re_proposed = log_re(proposed_b_mat, L, sds);
  // calculate the numerator
  vec numerator_b =
    logLik_long_proposed + logLik_surv_proposed + logLik_re_proposed;
  // log_rati
  vec log_ratio = numerator_b - denominator_b;
  uword n = log_ratio.n_elem;
  for (uword i = 0; i < n; i++) {
    if (std::isfinite(log_ratio.at(i)) && exp(log_ratio.at(i)) > R::runif(0, 1)) {
      acceptance_b.at(it, i) = 1;
      b_mat.row(i) = proposed_b_mat.row(i);
      logLik_long.at(i) = logLik_long_proposed.at(i);
      logLik_surv.at(i) = logLik_surv_proposed.at(i);
      logLik_re.at(i) = logLik_re_proposed.at(i);
      uword n_outcomes = eta.n_elem;
      for (uword j = 0; j < n_outcomes; j++) {
        eta.at(j).elem(idL_ind.at(i)) = eta_proposed.at(j).elem(idL_ind.at(i));
      }

      Wlong_H.rows(rows_Wlong_H.at(i)) =
        Wlong_H_proposed.rows(rows_Wlong_H.at(i));
      Wlong_h.row(i) = Wlong_h_proposed.row(i);
      Wlong_H2.rows(rows_Wlong_H.at(i)) =
        Wlong_H2_proposed.rows(rows_Wlong_H.at(i));

     // WlongH_alphas.rows(rows_Wlong_H.at(i)) =
      //  WlongH_alphas_proposed.rows(rows_Wlong_H.at(i));
     // Wlongh_alphas.row(i) = Wlongh_alphas_proposed.row(i);
     // WlongH2_alphas.rows(rows_Wlong_H.at(i)) =
     //   WlongH2_alphas_proposed.rows(rows_Wlong_H.at(i));
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

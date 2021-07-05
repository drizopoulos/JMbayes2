#ifndef JMBAYES2RE_H
#define JMBAYES2RE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_b (field<mat> &b, mat &b_mat, field<vec> &eta,
               vec &logLik_long, vec &logLik_surv, vec &logLik_re,
               mat &Wlong_H, mat &Wlong_h, mat &Wlong_H2,
               vec &WlongH_alphas, vec &Wlongh_alphas, vec &WlongH2_alphas,
               mat &scale_b,
               const field<uvec> &ind_RE,
               const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
               const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
               const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
               const mat &Wlong_bar, const mat &Wlong_sds, const field<vec> &betas,
               const vec &alphas, const uvec &id_H_, const uvec &id_h,
               const field<uvec> &FunForms, const field<uvec> &FunForms_ind,
               const List Funs_FunForms,
               const field<mat> &X, const field<mat> &Z,
               const field<uvec> &idL, const field<mat> &y,  const vec &sigmas,
               const vec &extra_parms, const CharacterVector &families,
               const CharacterVector &links, const field<uvec> &ids,
               const field<uvec> &unq_ids, const vec &W0H_bs_gammas,
               const vec &W0h_bs_gammas,
               const vec &W0H2_bs_gammas, const vec &WH_gammas,
               const vec &Wh_gammas, const vec &WH2_gammas,
               const vec &log_Pwk, const vec &log_Pwk2,
               const uvec &indFast_H, const uvec &indFast_h, const uvec &which_event,
               const uvec &which_right_event, const uvec &which_left,
               const uvec &which_interval, const bool &any_event,
               const bool &any_interval, const umat &ni_event,
               const mat &L, const vec &sds,
               const uword &it, mat &acceptance_b, cube &res_b, cube &res_b_last, 
               const bool &save_random_effects,
               const uword &n_burnin, const uword &n_iter, 
               const uword &GK_k, mat &cumsum_b, cube &outprod_b, 
               //
               const bool &recurrent, //!! new
               const vec &frailty_H, const vec &frailty_h, //!! new 
               const vec &alphaF_H, const vec &alphaF_h //!! new
) {
  uword n = b_mat.n_rows;
  uword nRE = b_mat.n_cols;
  // calculate denominator_b
  vec denominator_b = logLik_long + logLik_surv + logLik_re;
  for (uword j = 0; j < nRE; ++j) {
    // propose new random effects in mat and field<mat> form
    mat proposed_b_mat = propose_rnorm_mat(b_mat, scale_b, j);
    field<mat> proposed_b = mat2field(proposed_b_mat, ind_RE);

    // calculate log_lik_long based on proposed_b_mat
    field<vec> eta_proposed = linpred_mixed(X, betas, Z, proposed_b, idL);
    vec logLik_long_proposed = log_long(y, eta_proposed, sigmas, extra_parms,
                                        families, links, ids, unq_ids, n);
    // calculate Wlong_H, Wlong_h and Wlong_H2 based on the proposed_b
    // and calculate Wlong * alphas
    mat Wlong_H_proposed =
      calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, Wlong_sds, betas, proposed_b,
                      id_H_, FunForms, FunForms_ind, Funs_FunForms);
    vec WlongH_alphas_proposed = Wlong_H_proposed * alphas;

    mat Wlong_h_proposed(Wlong_h.n_rows, Wlong_h.n_cols);
    vec Wlongh_alphas_proposed(Wlongh_alphas.n_rows);
    if (any_event) {
      Wlong_h_proposed =
        calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, Wlong_sds, betas, proposed_b,
                        id_h, FunForms, FunForms_ind, Funs_FunForms);
      Wlongh_alphas_proposed = Wlong_h_proposed * alphas;
    }
    mat Wlong_H2_proposed(Wlong_H2.n_rows, Wlong_H2.n_cols);
    vec WlongH2_alphas_proposed(WlongH2_alphas.n_rows);
    if (any_interval) {
      Wlong_H2_proposed =
        calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, Wlong_sds, betas,
                        proposed_b, id_H_, FunForms, FunForms_ind, Funs_FunForms);
      WlongH2_alphas_proposed = Wlong_H2_proposed * alphas;
    }
    // calculate logLik_Surv_proposed
    vec logLik_surv_proposed =
      log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
               WH_gammas, Wh_gammas, WH2_gammas,
               WlongH_alphas_proposed, Wlongh_alphas_proposed,
               WlongH2_alphas_proposed,
               log_Pwk, log_Pwk2, indFast_H, indFast_h,
               which_event, which_right_event, which_left,
               any_interval, which_interval,
               recurrent, frailty_H, frailty_h, alphaF_H, alphaF_h); //!! new
    // logLik_re
    vec logLik_re_proposed = log_re(proposed_b_mat, L, sds);
    // calculate the numerator
    vec numerator_b =
      logLik_long_proposed + logLik_surv_proposed + logLik_re_proposed;
    // log_ratio
    vec log_ratio = numerator_b - denominator_b;
    for (uword i = 0; i < n; ++i) {
      double acc_i(0.0);
      if (std::isfinite(log_ratio.at(i)) &&
          exp(log_ratio.at(i)) > R::runif(0.0, 1.0)) {
        acc_i = 1.0;
        if (it > n_burnin - 1) acceptance_b.at(i, j) += 1.0;
        b_mat.row(i) = proposed_b_mat.row(i);
        denominator_b.at(i) = numerator_b.at(i);
        logLik_long.at(i) = logLik_long_proposed.at(i);
        logLik_surv.at(i) = logLik_surv_proposed.at(i);
        logLik_re.at(i) = logLik_re_proposed.at(i);
        uword first_H = GK_k * ni_event.at(i, 0);
        uword last_H = GK_k * ni_event.at(i, 1) - 1;
        Wlong_H.rows(first_H, last_H) = Wlong_H_proposed.rows(first_H, last_H);
        WlongH_alphas.rows(first_H, last_H) =
          WlongH_alphas_proposed.rows(first_H, last_H);
        if (any_event) {
          uword fitst_h = ni_event.at(i, 0);
          uword last_h = ni_event.at(i, 1) - 1;
          Wlong_h.rows(fitst_h, last_h) = Wlong_h_proposed.rows(fitst_h, last_h);
          Wlongh_alphas.rows(fitst_h, last_h) =
            Wlongh_alphas_proposed.rows(fitst_h, last_h);
        }
        if (any_interval) {
          Wlong_H2.rows(first_H, last_H) = Wlong_H2_proposed.rows(first_H, last_H);
          WlongH2_alphas.rows(first_H, last_H) =
            WlongH2_alphas_proposed.rows(first_H, last_H);
        }
      }
      if (it > 119) {
        scale_b.at(i, j) =
          robbins_monro(scale_b.at(i, j), acc_i, it - 100);
      }
    }
  }
  if (save_random_effects) {
    res_b.slice(it) = b_mat;
  } else {
    if (it > n_burnin - 1) {
      cumsum_b += b_mat;
      for (uword j = 0; j < n; j++) {
        outprod_b.slice(j) += b_mat.row(j).t() * b_mat.row(j);
      }
    }
  }
  if (it == n_iter - 1) {
    res_b_last.slice(0) = b_mat;
  }
  b = mat2field(b_mat, ind_RE);
  eta = linpred_mixed(X, betas, Z, b, idL);
}

#endif

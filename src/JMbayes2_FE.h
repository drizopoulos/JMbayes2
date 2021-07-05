#ifndef JMBAYES2FE_H
#define JMBAYES2FE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_betas (field<vec> &betas, mat &res_betas, field<vec> &acceptance_betas,
                   field<vec> &scale_betas, field<vec> &eta, vec &logLik_long,
                   vec &logLik_surv, mat &Wlong_H, mat &Wlong_h, mat &Wlong_H2,
                   vec &WlongH_alphas, vec &Wlongh_alphas, vec &WlongH2_alphas,
                   const vec &Tau_mean_betas_HC, const mat &prior_Tau_betas_HC,
                   const mat &b_mat, const mat &L, const vec &sds, const mat &X_dot,
                   const field<uvec> &ind_FE, // indices for the FE in res_betas[it,] belonging to the field betas. E.g., {{1,2,3}, {4, 5}, {6}}
                   const uvec &ind_FE_HC, // indices for the FE present in the HC (cols in res_betas)
                   const uvec &id_patt, // vector with the ids' outcome missing pattern
                   const field<uvec> &ind_RE_patt, // indices for the RE present in each outcome missing pattern (cols in D)
                   const field<uvec> &ind_FE_patt, // indices for the FE (in HC) present in each outcome missi
                   const uword &it,
                   const uvec &has_tilde_betas,
                   const field<mat> &X,
                   const field<mat> &Z,
                   const field<mat> &b,
                   const field<uvec> &idL,
                   const field<mat> &y,
                   const vec &sigmas,
                   const vec &extra_parms,
                   const CharacterVector &families,
                   const CharacterVector &links,
                   const field<uvec> &idL_lp_fast,
                   const field<vec> &prior_mean_betas_nHC,
                   field<mat> &prior_Tau_betas_nHC,
                   const field<uvec> &x_notin_z,
                   const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
                   const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
                   const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
                   const mat &Wlong_bar, const mat &Wlong_sds,
                   const uvec &id_H_, const uvec &id_h,
                   const field<uvec> &FunForms,
                   const field<uvec> &FunForms_ind,
                   const List Funs_FunForms,
                   const vec &alphas,
                   const bool &any_event, const bool &any_interval,
                   const vec &W0H_bs_gammas, const vec &W0h_bs_gammas, const vec &W0H2_bs_gammas,
                   const vec &WH_gammas, const vec &Wh_gammas, const vec &WH2_gammas,
                   const vec &log_Pwk, const vec &log_Pwk2,
                   const uvec &id_H_fast, const uvec &id_h_fast,
                   const uvec &which_event, const uvec &which_right_event, const uvec &which_left,
                   const uvec &which_interval, const field<uvec> &unq_idL,
                   const uword &n_burnin,
                   //
                   const bool &recurrent, //!! new
                   const vec &frailty_H, const vec &frailty_h, //!! new 
                   const vec &alphaF_H, const vec &alphaF_h //!! new
                  ) {
  uword n_b = b_mat.n_rows;
  // FE in HC - Gibbs sampling
  vec betas_vec = docall_rbindF(betas);
  uword patt_count = ind_RE_patt.n_elem; // number of unique outcome-missing patterns
  uword p_HC = ind_FE_HC.n_elem; // number of HC-FE
  uword q = b_mat.n_cols; // number of RE
  mat sum_JXDXJ(p_HC, p_HC, fill::zeros); // sum for the posterior's parameters
  vec sum_JXDu(p_HC, fill::zeros); // sum for the posterior's parameters
  mat U = L.each_row() % sds.t(); // RE vcov matrix Cholesky factorization (upper)
  field<mat> D_inv(patt_count); // all unique vcov_inv matrices across the missing outcome patterns
  for (uword i = 0; i < n_b; ++i) { // i-th patient
    /* I'm assuming that n_b in the total number of patients, including those
     * who only have survival outcome. I.e., b_mat has rows for them too.
     */
    if (i < patt_count && !ind_RE_patt.at(i).is_empty()) {
      /* obtain all unique vcov_inv matrices required for the sums in the posterior parameters
       * & jumps the pattern in which the patient misses all longitudinal outcomes
       */
      mat U_patt_inv = inv(trimatu(chol_update(U, ind_RE_patt.at(i))));
      D_inv.at(i) = U_patt_inv * U_patt_inv.t();
    }
    uword patt_i = id_patt.at(i); // id missing outcome pattern
    if (ind_FE_patt.at(patt_i).is_empty()) continue; // skip ids without longitudinal outcomes
    uvec ind_FE_i = ind_FE_patt.at(patt_i);
    uvec ind_RE_i = ind_RE_patt.at(patt_i);
    mat X_dot_i = X_dot.rows(i * q, (i + 1) * q - 1);
    X_dot_i = X_dot_i.submat(ind_RE_i, ind_FE_i);
    vec b_i = b_mat.row(i).t();
    vec u_i = b_i.rows(ind_RE_i) + X_dot_i * betas_vec.rows(ind_FE_i);
    mat D_inv_i = D_inv.at(patt_i);
    mat XD_i = X_dot_i.t() * D_inv_i;
    vec XDu = XD_i * u_i;
    mat XDX_i = XD_i * X_dot_i;
    sum_JXDu += add_zero_rows(XDu, p_HC, ind_FE_i);
    sum_JXDXJ += add_zero_colrows(XDX_i, p_HC, p_HC, ind_FE_i, ind_FE_i);
  }
  mat Sigma_1 = inv(prior_Tau_betas_HC + sum_JXDXJ); // improve via Cholesky decomposition
  vec mean_1 = Sigma_1 * (Tau_mean_betas_HC + sum_JXDu);
  mat U_1 = chol(Sigma_1);
  betas_vec.rows(ind_FE_HC) = propose_mvnorm_vec(U_1, 1.0) + mean_1;
  betas = vec2field(betas_vec, ind_FE);
  // update eta
  eta = linpred_mixed(X, betas, Z, b, idL);
  // update logLik_surv
  Wlong_H =
    calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, Wlong_sds, betas, b, id_H_,
                    FunForms, FunForms_ind, Funs_FunForms);
  WlongH_alphas = Wlong_H * alphas;
  if (any_event) {
    Wlong_h =
      calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, Wlong_sds, betas, b, id_h,
                      FunForms, FunForms_ind, Funs_FunForms);
    Wlongh_alphas = Wlong_h * alphas;
  }
  if (any_interval) {
    Wlong_H2 =
      calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, Wlong_sds, betas, b,
                      id_H_, FunForms, FunForms_ind, Funs_FunForms);
    WlongH2_alphas = Wlong_H2 * alphas;
  }
  logLik_surv = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                         WH_gammas, Wh_gammas, WH2_gammas,
                         WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                         log_Pwk, log_Pwk2, id_H_fast, id_h_fast,
                         which_event, which_right_event, which_left,
                         any_interval, which_interval,
                         recurrent, frailty_H, frailty_h, alphaF_H, alphaF_h); //!! new

  ///////////////////////////////////////////////////////////////////////////////
  // FE outside HC - Metropolis-Hastings sampling
  if (any(has_tilde_betas)) {
    uword n_outcomes = betas.n_elem;
    for (uword j = 0; j < n_outcomes; ++j) { // j-th outcome
      if (!has_tilde_betas.at(j)) continue; // skip outcomes without nHC-FE
      uvec ind_j = x_notin_z.at(j);
      uword n_betas = ind_j.n_rows;
      // denominator
      double sum_logLik_long_j =
        sum(log_long_i(y.at(j), eta.at(j), sigmas.at(j), extra_parms.at(j),
                       std::string(families[j]), std::string(links[j]),
                       idL_lp_fast.at(j)));
      vec ll(n_betas);
      /* improve: have input vec logPrior, and then use logPrior.at(j),
       *  if we accept logPrior.at(j) = logPrior_j_prop. To avoid re-calculations
       *  at each iteration
       */
      double logPrior_j =
        logPrior(betas.at(j).rows(ind_j), prior_mean_betas_nHC.at(j),
                 prior_Tau_betas_nHC.at(j), ll, 1.0, false);
      double denominator_j = sum_logLik_long_j + sum(logLik_surv) + logPrior_j;
      for (uword i = 0; i < n_betas; ++i) {
        // proposal
        field<vec> betas_prop = betas;
        betas_prop.at(j).rows(ind_j) =
          propose_norm(betas.at(j).rows(ind_j), scale_betas.at(j), i);
        double logPrior_j_prop =
          logPrior(betas_prop.at(j).rows(ind_j), prior_mean_betas_nHC.at(j),
                   prior_Tau_betas_nHC.at(j), ll, 1.0, false);
        // logLik_long proposal
        field<vec> eta_prop = linpred_mixed_i(eta, X, betas_prop, Z, b, idL, j);
        double sum_logLik_long_j_prop =
          sum(log_long_i(y.at(j), eta_prop.at(j), sigmas.at(j), extra_parms.at(j),
                         std::string(families[j]), std::string(links[j]),
                         idL_lp_fast.at(j)));
        // logLik_surv proposal
        mat Wlong_H_prop =
          calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, Wlong_sds, betas_prop, b,
                          id_H_, FunForms, FunForms_ind, Funs_FunForms);
        vec WlongH_alphas_prop = Wlong_H_prop * alphas;
        mat Wlong_h_prop(Wlong_h.n_rows, Wlong_h.n_cols);
        vec Wlongh_alphas_prop(Wlongh_alphas.n_rows);
        if (any_event) {
          Wlong_h_prop =
            calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, Wlong_sds, betas_prop,
                            b, id_h, FunForms, FunForms_ind, Funs_FunForms);
          Wlongh_alphas_prop = Wlong_h_prop * alphas;
        }
        mat Wlong_H2_prop(Wlong_H2.n_rows, Wlong_H2.n_cols);
        vec WlongH2_alphas_prop(WlongH2_alphas.n_rows);
        if (any_interval) {
          Wlong_H2_prop =
            calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, Wlong_sds, betas_prop,
                            b, id_H_, FunForms, FunForms_ind, Funs_FunForms);
          WlongH2_alphas_prop = Wlong_H2_prop * alphas;
        }
        vec logLik_surv_prop =
          log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                   WH_gammas, Wh_gammas, WH2_gammas,
                   WlongH_alphas_prop, Wlongh_alphas_prop, WlongH2_alphas_prop,
                   log_Pwk, log_Pwk2, id_H_fast, id_h_fast,
                   which_event, which_right_event, which_left,
                   any_interval, which_interval,
                   recurrent, frailty_H, frailty_h, alphaF_H, alphaF_h); //!! new
        // numerator
        double numerator_j =
          sum_logLik_long_j_prop + sum(logLik_surv_prop) + logPrior_j_prop;
        // Hastings ratio
        double acc_i(0.0);
        double log_ratio_j = numerator_j - denominator_j;
        if (std::isfinite(log_ratio_j) &&
            std::exp(log_ratio_j) > R::runif(0.0, 1.0)) {
          acc_i = 1.0;
          if (it > n_burnin - 1) acceptance_betas.at(j).at(i) += 1.0;
          betas.at(j) = betas_prop.at(j);
          eta = eta_prop;
          Wlong_H = Wlong_H_prop;
          Wlong_h = Wlong_h_prop;
          Wlong_H2 = Wlong_H2_prop;
          WlongH_alphas = WlongH_alphas_prop;
          Wlongh_alphas = Wlongh_alphas_prop;
          WlongH2_alphas = WlongH2_alphas_prop;
          logLik_surv = logLik_surv_prop;
          denominator_j = numerator_j;
        }
        if (it > 119) {
          scale_betas.at(j).at(i) =
            robbins_monro(scale_betas.at(j).at(i), acc_i, it - 100);
        }
      }
    }
  }
  // update logLik_long with all betas
  logLik_long = log_long(y, eta, sigmas, extra_parms, families, links,
                         idL_lp_fast, unq_idL, n_b);
  // save all results
  res_betas.row(it) = docall_rbindF(betas).t();
}

#endif

#ifndef JMBAYES2LOG_DENS
#define JMBAYES2LOG_DENS

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

double log_density_surv (const vec &W0H_bs_gammas,
                         const vec &W0h_bs_gammas,
                         const vec &W0H2_bs_gammas,
                         const vec &WH_gammas,
                         const vec &Wh_gammas,
                         const vec &WH2_gammas,
                         const vec &WlongH_alphas,
                         const vec &Wlongh_alphas,
                         const vec &WlongH2_alphas,
                         const vec &log_Pwk, const vec &log_Pwk2,
                         const uvec &indFast_H,
                         const uvec &which_event,
                         const uvec &which_right_event,
                         const uvec &which_left,
                         const bool &any_interval,
                         const uvec &which_interval) {
  vec lambda_H = W0H_bs_gammas + WH_gammas + WlongH_alphas;
  vec H = group_sum(exp(log_Pwk + lambda_H), indFast_H);
  int n = H.n_rows;
  vec lambda_h(n);
  lambda_h.elem(which_event) = W0h_bs_gammas.elem(which_event) +
    Wh_gammas.elem(which_event) + Wlongh_alphas.elem(which_event);
  vec log_Lik_surv(n);
  log_Lik_surv.elem(which_right_event) = - H.elem(which_right_event);
  log_Lik_surv.elem(which_event) += lambda_h.elem(which_event);
  log_Lik_surv.elem(which_left) = log1p(- exp(- H.elem(which_left)));
  vec lambda_H2(lambda_H.n_rows);
  vec H2(n);
  if (any_interval) {
    lambda_H2 = W0H2_bs_gammas + WH2_gammas + WlongH2_alphas;
    H2 = group_sum(exp(log_Pwk2 + lambda_H2), indFast_H);
    log_Lik_surv.elem(which_interval) = - H.elem(which_interval) +
      log(- expm1(- H2.elem(which_interval)));
  }
  double logLik = sum(log_Lik_surv);
  return logLik;
}

double logPC_D_sds (const vec &sds, const mat &L, const mat &b,
                    const double &prior_D_sds_df,
                    const double &prior_D_sds_sigma) {
  mat chol_Sigma = L.each_row() % sds.t();
  double log_p_b = sum(log_dmvnrm_chol(b, chol_Sigma));
  double log_p_sds = sum(log_dht(sds, prior_D_sds_sigma, prior_D_sds_df));
  double out = log_p_b + log_p_sds;
  return out;
}

double logPC_D_L (const mat &L, const vec &sds, const mat &b,
                  const double &prior_D_L_etaLKJ) {
  uword p = L.n_rows;
  mat chol_Sigma = L.each_row() % sds.t(); // check this
  double log_p_b = sum(log_dmvnrm_chol(b, chol_Sigma));
  double log_p_L(0.0);
  for (unsigned i = 1; i < p; ++i) {
    log_p_L += (p - i - 1.0 + 2.0 * prior_D_L_etaLKJ - 2.0) * log(L.at(i, i));
  }
  double out = log_p_b + log_p_L;
  return out;
}

#endif

#ifndef JMBAYES2LOGDENS_H
#define JMBAYES2LOGDENS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

vec log_long (const field<mat> &y, const field<vec> &eta, const vec &sigmas,
              const vec &extra_parms, const CharacterVector &families,
              const CharacterVector &links, const field<uvec> &ids,
              const field<uvec> &unq_ids) {
  uword n_outcomes = y.size();
  uvec ns(n_outcomes);
  for (uword i = 0; i < n_outcomes; ++i) {
    ns.at(i) = ids.at(i).n_rows;
  }
  uword n = ns.max();
  vec out(n, fill::zeros);
  for (uword i = 0; i < n_outcomes; ++i) {
    uvec id_i = ids.at(i);
    uvec unq_id_i = unq_ids.at(i);
    mat y_i = y.at(i);
    uword N = y_i.n_rows;
    vec log_contr(N);
    vec mu_i = mu_fun(eta.at(i), std::string(links[i]));
    double sigma_i = sigmas.at(i);
    double extr_prm_i = extra_parms.at(i);
    if (families[i] == "gaussian") {
      log_contr = log_dnorm(y_i, mu_i, sigma_i);
    } else if (families[i] == "Student-t") {
      log_contr = log_dt((y_i - mu_i) / sigma_i, extr_prm_i) - log(sigma_i);
    } else if (families[i] == "beta") {
      log_contr = log_dbeta(y_i, mu_i * sigma_i, sigma_i * (1 - mu_i));
    } else if (families[i] == "Gamma") {
      log_contr = log_dgamma(y_i, square(mu_i) / sigma_i, sigma_i / mu_i);
    } else if (families[i] == "unit Lindley") {
      vec theta = 1 / mu_i - 1;
      vec comp1 = 2 * log(theta) - log(1 + theta);
      vec comp2 = - 3 * log(1 - y_i);
      vec comp3 = - (theta * y_i) / (1 - y_i);
      log_contr = comp1 + comp2 + comp3;
    } else if (families[i] == "binomial") {
      uword k = y_i.n_cols;
      if (k == 2) {
        // y_i.col(1) has been set to y_i.col(0) + y_i.col(1)
        // in jm_fit(), i.e., y_i.col(1) is already the number of trials
        // not the number of failures
        log_contr = log_dbinom(y_i.col(0), y_i.col(1), mu_i);
      } else {
        log_contr = y_i % log(mu_i) + (1 - y_i) % log(1 - mu_i);
      }
    } else if (families[i] == "poisson") {
      log_contr = log_dpois(y_i, mu_i);
    } else if (families[i] == "negative binomial") {
      log_contr = log_dnbinom(y_i, mu_i, sigma_i);
    }  else if (families[i] == "beta binomial") {
      uword k = y_i.n_cols;
      if (k == 2) {
        // y_i.col(1) has been set to y_i.col(0) + y_i.col(1)
        // in jm_fit(), i.e., y_i.col(1) is already the number of trials
        // not the number of failures
        log_contr = log_dbbinom(y_i.col(0), y_i.col(1), mu_i, sigma_i);
      } else {
        vec ones(n, fill::ones);
        log_contr = log_dbbinom(y_i, ones, mu_i, sigma_i);
      }
    }
    out.elem(unq_id_i) += group_sum(log_contr, id_i);
  }
  return out;
}

vec log_surv (const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
                 const vec &W0H2_bs_gammas, const vec &WH_gammas,
                 const vec &Wh_gammas, const vec &WH2_gammas,
                 const vec &WlongH_alphas, const vec &Wlongh_alphas,
                 const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
                 const uvec &indFast_H, const uvec &which_event,
                 const uvec &which_right_event, const uvec &which_left,
                 const bool &any_interval, const uvec &which_interval) {
  vec lambda_H = W0H_bs_gammas + WH_gammas + WlongH_alphas;
  vec H = group_sum(exp(log_Pwk + lambda_H), indFast_H);
  uword n = H.n_rows;
  vec lambda_h(n);
  lambda_h.elem(which_event) = W0h_bs_gammas.elem(which_event) +
    Wh_gammas.elem(which_event) + Wlongh_alphas.elem(which_event);
  vec out(n);
  out.elem(which_right_event) = - H.elem(which_right_event);
  out.elem(which_event) += lambda_h.elem(which_event);
  out.elem(which_left) = log1p(- exp(- H.elem(which_left)));
  vec lambda_H2(lambda_H.n_rows);
  vec H2(n);
  if (any_interval) {
    lambda_H2 = W0H2_bs_gammas + WH2_gammas + WlongH2_alphas;
    H2 = group_sum(exp(log_Pwk2 + lambda_H2), indFast_H);
    out.elem(which_interval) = - H.elem(which_interval) +
      log(- expm1(- H2.elem(which_interval)));
  }
  return out;
}

vec log_re (const mat &b, const mat &L, const vec &sds) {
  mat chol_Sigma = L.each_row() % sds.t();
  vec out = log_dmvnrm_chol(b, chol_Sigma);
  return out;
}

vec logLik (const field<mat> &y, const field<vec> &eta, const vec &sigmas,
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
               const mat &b_mat, const mat &L, const vec &sds) {
  vec logLik_long = log_long(y, eta, sigmas, extra_parms, families, links, ids,
                             unq_ids);
  vec logLik_surv = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                             WH_gammas, Wh_gammas, WH2_gammas, WlongH_alphas,
                             Wlongh_alphas, WlongH2_alphas, log_Pwk,
                             log_Pwk2, indFast_H, which_event,
                             which_right_event, which_left, any_interval,
                             which_interval);
  vec logLik_re = log_re(b_mat, L, sds);
  vec out = logLik_long + logLik_surv + logLik_re;
  return out;
}

double logLik_prior (const mat &L, const vec &sds,
                     const double &prior_D_sds_df, const double &prior_D_sds_sigma,
                     const double &prior_D_L_etaLKJ,
                     const vec &bs_gammas, const vec &gammas, const vec &alphas,
                     const vec &prior_mean_bs_gammas, const mat &prior_Tau_bs_gammas,
                     const vec &prior_mean_gammas, const mat &prior_Tau_gammas,
                     const vec &prior_mean_alphas, const mat &prior_Tau_alphas,
                     const double &tau_bs_gammas,
                     double prior_A_tau_bs_gammas, double prior_B_tau_bs_gammas) {
  double out(0.0);
  out += sum(log_dht(sds, prior_D_sds_sigma, prior_D_sds_df));
  uword p = L.n_rows;
  double log_p_L(0.0);
  for (unsigned i = 1; i < p; ++i) {
    log_p_L += (p - i - 1.0 + 2.0 * prior_D_L_etaLKJ - 2.0) * log(L.at(i, i));
  }
  out += log_p_L;
  out += logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas,
                  tau_bs_gammas);
  out += logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0);
  out += logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);

  return out;
}

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

#endif

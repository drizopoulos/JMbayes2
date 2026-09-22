#ifndef JMBAYES2LOGDENS_H
#define JMBAYES2LOGDENS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

inline void log_long_i (const mat &y_i, const vec &eta_i, vec &mu_i,
                        const double sigma_i, const double &extr_prm_i,
                        const std::string &fam_i, const std::string &link_i,
                        const uvec &idFast_i, vec &log_contr, vec &out) {
    uword N = y_i.n_rows;
    mu_i = eta_i;
    mu_fun(mu_i, link_i);
    if (fam_i == "gaussian") {
        log_dnorm_void(y_i, mu_i, sigma_i, log_contr);
    } else if (fam_i == "Student's-t") {
        log_dt_void(y_i, mu_i, sigma_i, extr_prm_i, log_contr);
    } else if (fam_i == "beta") {
        log_dbeta_void(y_i, mu_i * sigma_i, sigma_i * (1.0 - mu_i), log_contr);
    } else if (fam_i == "Gamma") {
        log_dgamma_void(y_i, sigma_i, mu_i / sigma_i, log_contr);
    } else if (fam_i == "unit Lindley") {
        const double* yy = y_i.memptr();
        const double* mu_ptr = mu_i.memptr();
        double* out_ptr = log_contr.memptr();
        for (uword i = 0; i < N; ++i) {
            double theta = 1.0 / mu_ptr[i] - 1.0;
            double y_val = yy[i];
            out_ptr[i] = 2.0 * std::log(theta) - std::log(1.0 + theta)
                - 3.0 * std::log(1.0 - y_val)
                - (theta * y_val) / (1.0 - y_val);
        }
    } else if (fam_i == "censored normal") {
        const double* yy = y_i.colptr(0);
        const double* cens = y_i.colptr(1);
        const double* mu_ptr = mu_i.memptr();
        double* out_ptr = log_contr.memptr();
        for (uword i = 0; i < N; ++i) {
            if (cens[i] == 0.0) {
                out_ptr[i] = R::dnorm(yy[i], mu_ptr[i], sigma_i, 1);
            } else if (cens[i] == 1.0) {
                out_ptr[i] = R::pnorm(yy[i], mu_ptr[i], sigma_i, 1, 1);
            } else { // cens == 2.0
                out_ptr[i] = R::pnorm(yy[i], mu_ptr[i], sigma_i, 0, 1);
            }
        }
    } else if (fam_i == "binomial") {
        if (y_i.n_cols == 2) {
            log_dbinom_void(y_i.col(0), y_i.col(1), mu_i, log_contr);
        } else {
            log_dbernoulli_void(y_i, mu_i, log_contr);
        }
    } else if (fam_i == "poisson") {
        log_dpois_void(y_i, mu_i, log_contr);
    } else if (fam_i == "negative binomial") {
        log_dnbinom_void(y_i, mu_i, sigma_i, log_contr);
    } else if (fam_i == "beta binomial") {
        if (y_i.n_cols == 2) {
            log_dbbinom_void(y_i.col(0), y_i.col(1), mu_i, sigma_i, log_contr);
        } else {
            vec ones(N, fill::ones);
            log_dbbinom_void(y_i, ones, mu_i, sigma_i, log_contr);
        }
    }
    group_sum(log_contr, idFast_i, out);
}

inline void log_long (const field<mat> &y, const field<vec> &eta, field<vec> &mu,
                      const vec &sigmas, const vec &extra_parms,
                      const std::vector<std::string> &families,
                      const std::vector<std::string> &links,
                      const field<uvec> &idFast, const field<uvec> &unq_ids,
                      vec &out, field<vec> &log_contr_obs_workspace,
                      field<vec> &log_contr_subj_workspace) {
    uword n_outcomes = y.size();
    out.zeros();
    for (uword i = 0; i < n_outcomes; ++i) {
        const mat& y_i = y.at(i);
        const vec& eta_i = eta.at(i);
        double sigma_i = sigmas.at(i);
        double extr_prm_i = extra_parms.at(i);
        const std::string& fam_i = families[i];
        const std::string& link_i = links[i];
        const uvec& idFast_i = idFast.at(i);
        const uvec& unq_id_i = unq_ids.at(i);
        log_long_i(y_i, eta_i, mu.at(i), sigma_i, extr_prm_i, fam_i, link_i,
                   idFast_i, log_contr_obs_workspace.at(i),
                   log_contr_subj_workspace.at(i));
        out.rows(unq_id_i) += log_contr_subj_workspace.at(i);
    }
}

inline void log_surv (const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
              const vec &W0H2_bs_gammas, const vec &WH_gammas,
              const vec &Wh_gammas, const vec &WH2_gammas,
              const vec &WlongH_alphas, const vec &Wlongh_alphas,
              const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
              const vec &log_weights, const uvec &ind_h2, const uvec &intgr_ind, const bool &intgr,
              const uvec &indFast_H, const uvec &indFast_h, const uvec &which_event,
              const uvec &which_right_event, const uvec &which_left,
              const bool &any_interval, const uvec &which_interval,
              const bool &recurrent,
              const vec &frailtyH_sigmaF_alphaF,
              const vec &frailtyh_sigmaF_alphaF,
              vec &lambda_H, vec &H, vec &lambda_H2, vec &H2, vec &out,
              vec &logLik_surv) {
    lambda_H = log_Pwk + W0H_bs_gammas + WH_gammas + WlongH_alphas;
    if (recurrent) {
        lambda_H += frailtyH_sigmaF_alphaF;
    }
    lambda_H = exp(lambda_H);
    group_sum(lambda_H, indFast_H, H);
    out.rows(which_right_event) = -H.rows(which_right_event);
    if (which_event.n_elem > 0) {
        out.rows(which_event) += log_weights.rows(which_event) +
            W0h_bs_gammas.rows(which_event) +
            Wh_gammas.rows(which_event) +
            Wlongh_alphas.rows(which_event);
        if (recurrent) {
            out.rows(which_event) += frailtyh_sigmaF_alphaF;
        }
    }
    if (which_left.n_elem > 0) {
        out.rows(which_left) = log1p(-exp(-H.rows(which_left)));
    }
    if (any_interval) {
        lambda_H2 = log_Pwk2 + W0H2_bs_gammas + WH2_gammas + WlongH2_alphas;
        lambda_H2 = exp(lambda_H2);
        group_sum(lambda_H2, indFast_H, H2);
        out.rows(which_interval) = -H.rows(which_interval) +
            log(-expm1(-H2.rows(which_interval)));
    }
    if (intgr) {
        out = lse(out, ind_h2, intgr_ind);
    }
    group_sum(out, indFast_h, logLik_surv);
}

inline void log_surv2 (const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
                      const vec &W0H2_bs_gammas, const vec &WH_gammas,
                      const vec &Wh_gammas, const vec &WH2_gammas,
                      const vec &WlongH_alphas, const vec &Wlongh_alphas,
                      const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
                      const vec &log_weights, const uvec &ind_h2, const uvec &intgr_ind, const bool &intgr,
                      const uvec &indFast_H, const uvec &indFast_h, const uvec &which_event,
                      const uvec &which_right_event, const uvec &which_left,
                      const bool &any_interval, const uvec &which_interval,
                      const bool &recurrent,
                      const vec &frailtyH_sigmaF_alphaF,
                      const vec &frailtyh_sigmaF_alphaF,
                      vec &lambda_H, vec &H, vec &lambda_H2, vec &H2, vec &out,
                      vec &logLik_surv) {
    // 1. Calculate lambda_H (Fused Math)
    uword N_H = log_Pwk.n_elem;
    double* lam_H_ptr = lambda_H.memptr();
    const double* log_Pwk_ptr = log_Pwk.memptr();
    const double* W0H_ptr = W0H_bs_gammas.memptr();
    const double* WH_ptr = WH_gammas.memptr();
    const double* WlongH_ptr = WlongH_alphas.memptr();
    if (recurrent) {
        const double* frail_ptr = frailtyH_sigmaF_alphaF.memptr();
        for (uword r = 0; r < N_H; ++r) {
            lam_H_ptr[r] = std::exp(log_Pwk_ptr[r] + W0H_ptr[r] + WH_ptr[r] +
                WlongH_ptr[r] + frail_ptr[r]);
        }
    } else {
        for (uword r = 0; r < N_H; ++r) {
            lam_H_ptr[r] = std::exp(log_Pwk_ptr[r] + W0H_ptr[r] + WH_ptr[r] +
                WlongH_ptr[r]);
        }
    }
    group_sum(lambda_H, indFast_H, H);
    double* out_ptr = out.memptr();
    const double* H_ptr = H.memptr();
    // 2. Right Censored Events
    uword N_right = which_right_event.n_elem;
    const uword* right_ptr = which_right_event.memptr();
    for(uword i = 0; i < N_right; ++i) {
        uword idx = right_ptr[i];
        out_ptr[idx] = -H_ptr[idx];
    }
    // 3. Observed Events
    uword N_event = which_event.n_elem;
    if (N_event > 0) {
        const uword* event_ptr = which_event.memptr();
        const double* log_w_ptr = log_weights.memptr();
        const double* W0h_ptr = W0h_bs_gammas.memptr();
        const double* Wh_ptr = Wh_gammas.memptr();
        const double* Wlongh_ptr = Wlongh_alphas.memptr();
        if (recurrent) {
            const double* frail_h_ptr = frailtyh_sigmaF_alphaF.memptr();
            for(uword i = 0; i < N_event; ++i) {
                uword idx = event_ptr[i];
                out_ptr[idx] += log_w_ptr[idx] + W0h_ptr[idx] + Wh_ptr[idx] +
                    Wlongh_ptr[idx] + frail_h_ptr[idx];
            }
        } else {
            for(uword i = 0; i < N_event; ++i) {
                uword idx = event_ptr[i];
                out_ptr[idx] += log_w_ptr[idx] + W0h_ptr[idx] + Wh_ptr[idx] +
                    Wlongh_ptr[idx];
            }
        }
    }
    // 4. Left Censored Events
    uword N_left = which_left.n_elem;
    if (N_left > 0) {
        const uword* left_ptr = which_left.memptr();
        for(uword i = 0; i < N_left; ++i) {
            uword idx = left_ptr[i];
            out_ptr[idx] = std::log1p(-std::exp(-H_ptr[idx]));
        }
    }
    // 5. Interval Censored Events
    if (any_interval) {
        uword N_H2 = log_Pwk2.n_elem;
        double* lam_H2_ptr = lambda_H2.memptr();
        const double* log_Pwk2_ptr = log_Pwk2.memptr();
        const double* W0H2_ptr = W0H2_bs_gammas.memptr();
        const double* WH2_ptr = WH2_gammas.memptr();
        const double* WlongH2_ptr = WlongH2_alphas.memptr();
        for (uword r = 0; r < N_H2; ++r) {
            lam_H2_ptr[r] = std::exp(log_Pwk2_ptr[r] + W0H2_ptr[r] + WH2_ptr[r] + WlongH2_ptr[r]);
        }
        group_sum(lambda_H2, indFast_H, H2);
        uword N_interval = which_interval.n_elem;
        const uword* int_ptr = which_interval.memptr();
        const double* H2_ptr = H2.memptr();

        for(uword i = 0; i < N_interval; ++i) {
            uword idx = int_ptr[i];
            out_ptr[idx] = -H_ptr[idx] + std::log(-std::expm1(-H2_ptr[idx]));
        }
    }
    // 6. Final Integration and Aggregation
    if (intgr) {
        out = lse(out, ind_h2, intgr_ind);
    }

    group_sum(out, indFast_h, logLik_surv);
}

vec log_surv_old (const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
              const vec &W0H2_bs_gammas, const vec &WH_gammas,
              const vec &Wh_gammas, const vec &WH2_gammas,
              const vec &WlongH_alphas, const vec &Wlongh_alphas,
              const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
              const vec &log_weights, const uvec &ind_h2, const uvec &intgr_ind, const bool &intgr,
              const uvec &indFast_H, const uvec &indFast_h, const uvec &which_event,
              const uvec &which_right_event, const uvec &which_left,
              const bool &any_interval, const uvec &which_interval,
              const bool &calculate_sum = true) {
  vec lambda_H = W0H_bs_gammas + WH_gammas + WlongH_alphas;
  vec H = group_sum(exp(log_Pwk + lambda_H), indFast_H);
  uword n = H.n_rows;
  vec lambda_h = log_weights;
  lambda_h.rows(which_event) = W0h_bs_gammas.rows(which_event) +
    Wh_gammas.rows(which_event) + Wlongh_alphas.rows(which_event);
  vec out(n);
  out.rows(which_right_event) = - H.rows(which_right_event);
  out.rows(which_event) += lambda_h.rows(which_event);
  out.rows(which_left) = log1p(- exp(- H.rows(which_left)));
  vec lambda_H2(lambda_H.n_rows);
  vec H2(n);
  if (any_interval) {
    lambda_H2 = W0H2_bs_gammas + WH2_gammas + WlongH2_alphas;
    H2 = group_sum(exp(log_Pwk2 + lambda_H2), indFast_H);
    out.rows(which_interval) = - H.rows(which_interval) +
      log(- expm1(- H2.rows(which_interval)));
  }
  if (intgr) {
      out = lse(out, ind_h2, intgr_ind);
  }
  if (calculate_sum) {
      out = group_sum(out, indFast_h);
  }
  return out;
}

inline vec log_re (const mat &b, const mat &L, const vec &sds) {
    uword k = b.n_cols;
    double log_det = -arma::sum(arma::log(L.diag())) - arma::sum(arma::log(sds));
    double constants = -(double)k / 2.0 * log2pi;
    double other_terms = constants + log_det;
    mat B_scaled = b.each_row() / sds.t();
    mat Z_transposed = arma::solve(arma::trimatl(L.t()), B_scaled.t());
    vec sq_dist = arma::sum(arma::square(Z_transposed), 0).t();
    return other_terms - 0.5 * sq_dist;
}

inline void log_re_onlyRE (const mat &b_prop, const mat &V_Sigma,
                           double other_terms, mat &Z_workspace, vec &log_re) {
    // Fast Triangular Multiplication (TRMM)
    // V_Sigma is upper-triangular, so trimatu() safely bypasses half the math.
    Z_workspace = b_prop * arma::trimatu(V_Sigma);
    // Calculate squared distances for all subjects
    // sum(..., 1) sums across the rows, giving an (n x 1) vector
    log_re = other_terms - 0.5 * arma::sum(arma::square(Z_workspace), 1);
}

inline void log_re_onlySDS (const mat &b, const mat &V_R, double log_det_V_R,
                           const vec &sds, mat &B_scaled_workspace,
                           mat &Z_workspace, vec &log_re) {
    uword k = b.n_cols;
    // Calculate constants
    double log_det_V_Sigma = log_det_V_R - arma::sum(arma::log(sds));
    double other_terms = -(double)k / 2.0 * log2pi + log_det_V_Sigma;
     // Fast Triangular Multiplication (TRMM)
     B_scaled_workspace = b.each_row() / sds.t();
    // By wrapping V_R in trimatu(), we tell the BLAS backend to skip half the math.
    Z_workspace = B_scaled_workspace * arma::trimatu(V_R);
    // Final vectorized computation
    log_re = other_terms - 0.5 * arma::sum(arma::square(Z_workspace), 1);
}

inline void log_re_onlyL (const mat &b_scaled, const mat &L, double sum_log_sds,
                         mat &Z_transposed, vec &log_re) {
    uword k = b_scaled.n_cols;
    // Calculate constants
    double log_det = -arma::sum(arma::log(L.diag())) - sum_log_sds;
    double constants = -(double)k / 2.0 * log2pi;
    double other_terms = constants + log_det;
    // Fast Triangular Solve (TRSM)
    // We want Z = b_scaled * L^-1.
    // We achieve this safely by solving L^T * Z^T = b_scaled^T
    Z_transposed = arma::solve(arma::trimatl(L.t()), b_scaled.t());
    // Calculate squared distances for all subjects
    // sum(..., 0) sums down the columns, returning a (1 x n) row vector.
    // .t() flips it to an (n x 1) column vector to match your expected output.
    // Final vectorized computation
    log_re = other_terms - 0.5 * arma::sum(arma::square(Z_transposed), 0).t();
}

/*
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
 */

vec logLik_jm_stripped (
    const field<vec> &betas, const field<mat> &b, const vec &sigmas,
    const vec &bs_gammas, const vec &gammas, const vec &alphas, const vec &tau_bs_gammas,
    const mat &L, const vec &sds,
    /////////////
    const field<mat> &y, const field<mat> &X, const field<mat> &Xbar, const field<mat> &Z,
    const vec &extra_parms, const std::vector<std::string> &families,
    const std::vector<std::string> &links,
    const field<uvec> &idL, const field<uvec> &idL_lp_fast, const field<uvec> &unq_idL,
    /////////////
    const mat &W0_H, const mat &W0_h, const mat &W0_H2,
    const mat &W_H, const mat &W_h, const mat &W_H2,
    const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
    const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
    const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
    const mat &Wlong_bar, const mat &Wlong_sds, const mat &W_sds,
    const bool &any_event, const bool &any_interval, const bool &any_gammas,
    const field<uvec> &FunForms, const List &Funs_FunForms,
    const uvec &id_H_, const uvec &id_h,
    const vec &log_Pwk, const vec &log_Pwk2, const vec &log_weights,
    const uvec &id_h2, const uvec &intgr_ind, const bool &intgr,
    const uvec &id_H_fast, const uvec &id_h_fast,
    const uvec &which_event, const uvec &which_right_event,
    const uvec &which_left, const uvec &which_interval,
    const bool &recurrent, const vec &alphaF, const vec &frailty,
    const field<uvec> &which_term_H, const field<uvec> &which_term_h, const bool &any_terminal,
    const vec &sigmaF, vec &lambda_H_workspace, vec &H_workspace,
    vec &lambda_H2_workspace, vec &H2_workspace, vec &surv_out_workspace,
    vec &logLik_long, vec &logLik_surv,
    field<vec> &mu_obs_workspace,
    field<vec> &log_contr_obs_workspace,
    field<vec> &log_contr_subj_workspace) {
  //uword n = b.at(0).n_rows;
  /////////////
  field<vec> betas_ = betas;
  // set intercept to centered covariates
  for (uword j = 0; j < Xbar.n_elem; ++j) {
    betas_.at(j).at(0) += as_scalar(Xbar.at(j) * betas.at(j));
  }
  field<vec> eta = linpred_mixed(X, betas_, Z, b, idL);
  log_long(y, eta, mu_obs_workspace, sigmas, extra_parms, families, links,
           idL_lp_fast, unq_idL, logLik_long, log_contr_obs_workspace,
           log_contr_subj_workspace);
  /////////////
  vec W0H_bs_gammas = W0_H * bs_gammas;
  vec W0h_bs_gammas(W0_h.n_rows);
  vec W0H2_bs_gammas(W0_H2.n_rows);
  if (any_event) {
    W0h_bs_gammas = W0_h * bs_gammas;
  }
  if (any_interval) {
    W0H2_bs_gammas = W0_H2 * bs_gammas;
  }
  vec WH_gammas(W0_H.n_rows);
  vec Wh_gammas(W0_h.n_rows);
  vec WH2_gammas(W0_H2.n_rows);
  vec gammas_ = gammas % W_sds.t();
  if (any_gammas) {
    WH_gammas = W_H * gammas_;
  }
  if (any_gammas && any_event) {
    Wh_gammas = W_h * gammas_;
  }
  if (any_gammas && any_interval) {
    WH2_gammas = W_H2 * gammas_;
  }
  mat Wlong_H =
    calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, Wlong_sds, betas_, b, id_H_,
                    FunForms, Funs_FunForms);
  vec alphas_ = alphas % Wlong_sds.t();
  vec WlongH_alphas = Wlong_H * alphas_;
  mat Wlong_h(W0_h.n_rows, alphas.n_rows);
  vec Wlongh_alphas(W0_h.n_rows);
  if (any_event) {
    Wlong_h =
      calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, Wlong_sds, betas_, b, id_h,
                      FunForms, Funs_FunForms);
    Wlongh_alphas = Wlong_h * alphas_;
  }
  mat Wlong_H2(W0_H2.n_rows, alphas.n_rows);
  vec WlongH2_alphas(W0_H2.n_rows);
  if (any_interval) {
    Wlong_H2 =
      calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, Wlong_sds, betas_, b, id_H_,
                      FunForms, Funs_FunForms);
    WlongH2_alphas = Wlong_H2 * alphas_;
  }
  vec alphaF_H(WH_gammas.n_rows, fill::ones);
  vec alphaF_h(Wh_gammas.n_rows, fill::ones);
  if(any_terminal) {
    for (uword j = 0; j < alphaF.n_rows; ++j) {
      alphaF_H.rows(which_term_H.at(j)).fill(alphaF.at(j));
      alphaF_h.rows(which_term_h.at(j)).fill(alphaF.at(j));
    }
  }
  vec frailty_H(WH_gammas.n_rows, fill::zeros);
  vec frailty_h(Wh_gammas.n_rows, fill::zeros);
  frailty_h = frailty.rows(id_h);
  frailty_H = frailty.rows(id_H_);
  vec frailtyH_sigmaF_alphaF(WH_gammas.n_rows, fill::zeros);
  vec frailtyh_sigmaF_alphaF(which_event.n_rows, fill::zeros);
  frailtyH_sigmaF_alphaF = frailty_H % alphaF_H * sigmaF;
  frailtyh_sigmaF_alphaF = frailty_h.rows(which_event) % alphaF_h.rows(which_event) * sigmaF;
  log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
           WH_gammas, Wh_gammas, WH2_gammas,
           WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
           log_Pwk, log_Pwk2, log_weights, id_h2, intgr_ind, intgr,
           id_H_fast, id_h_fast,
           which_event, which_right_event, which_left,
           any_interval, which_interval,
           recurrent, frailtyH_sigmaF_alphaF, frailtyh_sigmaF_alphaF,
           lambda_H_workspace, H_workspace,
           lambda_H2_workspace, H2_workspace, surv_out_workspace,
           logLik_surv);
  mat b_mat = docall_cbindF(b);
  vec logLik_re = log_re(b_mat, L, sds);
  vec out = logLik_long + logLik_surv + logLik_re;
  if(recurrent) {
    vec logLik_frailty = log_dnorm(frailty, vec(frailty.n_elem, fill::zeros), 1.0);
    out += logLik_frailty;
  }
  return out;
}

#endif

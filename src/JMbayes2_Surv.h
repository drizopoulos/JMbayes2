#ifndef JMBAYES2SURV_H
#define JMBAYES2SURV_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

double logPrior_surv (
    const vec &bs_gammas, const vec&gammas, const vec &alphas,
    const field<vec> &prior_mean_bs_gammas, field<mat> &prior_Tau_bs_gammas,
    const vec &tau_bs_gammas,
    const vec &prior_mean_gammas, mat &prior_Tau_gammas, const vec &lambda_gammas,
    const double &tau_gammas, const bool &shrink_gammas,
    const vec &prior_mean_alphas, mat &prior_Tau_alphas, const vec &lambda_alphas,
    const double &tau_alphas, const bool &shrink_alphas,
    const bool &recurrent, const vec &alphaF, const vec prior_mean_alphaF,
    mat &prior_Tau_alphaF, const vec &lambda_alphaF, const double &tau_alphaF, 
    const bool &shrink_alphaF) {
  uword n_strata = prior_mean_bs_gammas.n_elem;
  uword n_per_stratum = bs_gammas.n_rows / n_strata;
  double out(0.0);
  for (uword i = 0; i < n_strata; ++i) {
    vec mu = prior_mean_bs_gammas.at(i);
    out += logPrior(bs_gammas.rows(i * n_per_stratum, (i + 1) * n_per_stratum - 1),
                    mu, prior_Tau_bs_gammas.at(i), mu.ones(), tau_bs_gammas.at(i),
                    false);
  }
  out += logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, lambda_gammas,
                  tau_gammas, shrink_gammas);
  out += logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, lambda_alphas,
                  tau_alphas, shrink_alphas);
  if (recurrent) {
    out += logPrior(alphaF, prior_mean_alphaF, prior_Tau_alphaF, lambda_alphaF,
                    tau_alphaF, shrink_alphaF);
  }
  return out;
}


void update_bs_gammas (vec &bs_gammas, const vec &gammas, const vec &alphas,
                       vec &W0H_bs_gammas, vec &W0h_bs_gammas, vec &W0H2_bs_gammas,
                       const vec &WH_gammas, const vec &Wh_gammas, const vec &WH2_gammas,
                       const vec &WlongH_alphas, const vec &Wlongh_alphas, const vec &WlongH2_alphas,
                       const vec &log_Pwk, const vec &log_Pwk2,
                       const uvec &indFast_H, const uvec &indFast_h,
                       const uvec &which_event, const uvec &which_right_event,
                       const uvec &which_left, const uvec &which_interval,
                       const bool &any_event, const bool &any_interval,
                       const field<vec> &prior_mean_bs_gammas, field<mat> &prior_Tau_bs_gammas,
                       const vec &tau_bs_gammas,
                       const vec &prior_mean_gammas, mat &prior_Tau_gammas,
                       const vec &lambda_gammas, const double &tau_gammas, const bool &shrink_gammas,
                       const vec &prior_mean_alphas, mat &prior_Tau_alphas,
                       const vec &lambda_alphas, const double &tau_alphas, const bool &shrink_alphas,
                       vec &logLik_surv, double &denominator_surv, const uword &it,
                       /////
                       const mat &W0_H, const mat &W0_h, const mat &W0_H2,
                       vec &scale_bs_gammas, mat &acceptance_bs_gammas,
                       mat &res_bs_gammas,
                       const bool &recurrent,
                       const vec &frailtyH_sigmaF_alphaF, const vec &frailtyh_sigmaF_alphaF,
                       const vec &alphaF, const vec prior_mean_alphaF,
                       mat &prior_Tau_alphaF, const vec &lambda_alphaF,
                       const double &tau_alphaF, const bool &shrink_alphaF) {
  for (uword i = 0; i < bs_gammas.n_rows; ++i) {
    vec proposed_bs_gammas = propose_norm(bs_gammas, scale_bs_gammas, i);
    vec proposed_W0H_bs_gammas = W0_H * proposed_bs_gammas;
    vec proposed_W0h_bs_gammas(W0_h.n_rows);
    vec proposed_W0H2_bs_gammas(W0_H2.n_rows);
    if (any_event) {
      proposed_W0h_bs_gammas = W0_h * proposed_bs_gammas;
    }
    if (any_interval) {
      proposed_W0H2_bs_gammas = W0_H2 * proposed_bs_gammas;
    }
    vec logLik_surv_proposed =
      log_surv(proposed_W0H_bs_gammas, proposed_W0h_bs_gammas, proposed_W0H2_bs_gammas,
               WH_gammas, Wh_gammas, WH2_gammas,
               WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
               log_Pwk, log_Pwk2, indFast_H, indFast_h,
               which_event, which_right_event, which_left,
               any_interval, which_interval,
               recurrent, frailtyH_sigmaF_alphaF, frailtyh_sigmaF_alphaF);
    double numerator_surv =
      sum(logLik_surv_proposed) +
      logPrior_surv(proposed_bs_gammas, gammas, alphas, prior_mean_bs_gammas,
                    prior_Tau_bs_gammas, tau_bs_gammas,
                    prior_mean_gammas, prior_Tau_gammas, lambda_gammas, tau_gammas, shrink_gammas,
                    prior_mean_alphas, prior_Tau_alphas, lambda_alphas, tau_alphas, shrink_alphas,
                    recurrent, alphaF, prior_mean_alphaF, prior_Tau_alphaF,
                    lambda_alphaF, tau_alphaF, shrink_alphaF);
    double log_ratio = numerator_surv - denominator_surv;
    if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0.0, 1.0)) {
      bs_gammas = proposed_bs_gammas;
      W0H_bs_gammas = proposed_W0H_bs_gammas;
      if (any_event) {
        W0h_bs_gammas = proposed_W0h_bs_gammas;
      }
      if (any_interval) {
        W0H2_bs_gammas = proposed_W0H2_bs_gammas;
      }
      logLik_surv = logLik_surv_proposed;
      denominator_surv = numerator_surv;
      acceptance_bs_gammas.at(it, i) = 1;
    }
    if (it > 19) {
      scale_bs_gammas.at(i) =
        robbins_monro(scale_bs_gammas.at(i),
                      acceptance_bs_gammas.at(it, i), it);
    }
    res_bs_gammas.at(it, i) = bs_gammas.at(i);
  }
}

void update_gammas (const vec &bs_gammas, vec &gammas, const vec &alphas,
                    const vec &W0H_bs_gammas, const vec &W0h_bs_gammas, const vec &W0H2_bs_gammas,
                    vec &WH_gammas, vec &Wh_gammas, vec &WH2_gammas,
                    const vec &WlongH_alphas, const vec &Wlongh_alphas, const vec &WlongH2_alphas,
                    const vec &log_Pwk, const vec &log_Pwk2,
                    const uvec &indFast_H, const uvec &indFast_h,
                    const uvec &which_event, const uvec &which_right_event,
                    const uvec &which_left, const uvec &which_interval,
                    const bool &any_event, const bool &any_interval,
                    const field<vec> &prior_mean_bs_gammas, field<mat> &prior_Tau_bs_gammas,
                    const vec &tau_bs_gammas,
                    const vec &prior_mean_gammas, mat &prior_Tau_gammas,
                    const vec &lambda_gammas, const double &tau_gammas, const bool &shrink_gammas,
                    const vec &prior_mean_alphas, mat &prior_Tau_alphas,
                    const vec &lambda_alphas, const double &tau_alphas, const bool &shrink_alphas,
                    vec &logLik_surv, double &denominator_surv, const uword &it,
                    /////
                    const mat &W_H, const mat &W_h, const mat &W_H2,
                    vec &scale_gammas, mat &acceptance_gammas, mat &res_gammas,
                    const bool &recurrent,
                    const vec &frailtyH_sigmaF_alphaF, const vec &frailtyh_sigmaF_alphaF,
                    const vec &alphaF, const vec prior_mean_alphaF,
                    mat &prior_Tau_alphaF, const vec &lambda_alphaF,
                    const double &tau_alphaF, const bool &shrink_alphaF) {
  for (uword i = 0; i < gammas.n_rows; ++i) {
    vec proposed_gammas = propose_norm(gammas, scale_gammas, i);
    vec proposed_WH_gammas = W_H * proposed_gammas;
    vec proposed_Wh_gammas(W_h.n_rows);
    vec proposed_WH2_gammas(W_H2.n_rows);
    if (any_event) {
      proposed_Wh_gammas = W_h * proposed_gammas;
    }
    if (any_interval) {
      proposed_WH2_gammas = W_H2 * proposed_gammas;
    }
    vec logLik_surv_proposed =
      log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
               proposed_WH_gammas, proposed_Wh_gammas, proposed_WH2_gammas,
               WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
               log_Pwk, log_Pwk2, indFast_H, indFast_h,
               which_event, which_right_event, which_left,
               any_interval, which_interval,
               recurrent, frailtyH_sigmaF_alphaF, frailtyh_sigmaF_alphaF);
    double numerator_surv =
      sum(logLik_surv_proposed) +
      logPrior_surv(bs_gammas, proposed_gammas, alphas, prior_mean_bs_gammas,
                    prior_Tau_bs_gammas, tau_bs_gammas,
                    prior_mean_gammas, prior_Tau_gammas, lambda_gammas, tau_gammas, shrink_gammas,
                    prior_mean_alphas, prior_Tau_alphas, lambda_alphas, tau_alphas, shrink_alphas,
                    recurrent, alphaF, prior_mean_alphaF, prior_Tau_alphaF,
                    lambda_alphaF, tau_alphaF, shrink_alphaF);
    double log_ratio = numerator_surv - denominator_surv;
    if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0.0, 1.0)) {
      gammas = proposed_gammas;
      WH_gammas = proposed_WH_gammas;
      if (any_event) {
        Wh_gammas = proposed_Wh_gammas;
      }
      if (any_interval) {
        WH2_gammas = proposed_WH2_gammas;
      }
      logLik_surv = logLik_surv_proposed;
      denominator_surv = numerator_surv;
      acceptance_gammas.at(it, i) = 1;
    }
    if (it > 19) {
      scale_gammas.at(i) =
        robbins_monro(scale_gammas.at(i),
                      acceptance_gammas.at(it, i), it);
    }
    // store results
    res_gammas.at(it, i) = gammas.at(i);
  }
}

void update_alphas (const vec &bs_gammas, const vec &gammas, vec &alphas,
                    const vec &W0H_bs_gammas, const vec &W0h_bs_gammas, const vec &W0H2_bs_gammas,
                    const vec &WH_gammas, const vec &Wh_gammas, const vec &WH2_gammas,
                    vec &WlongH_alphas, vec &Wlongh_alphas, vec &WlongH2_alphas,
                    const vec &log_Pwk, const vec &log_Pwk2,
                    const uvec &indFast_H, const uvec &indFast_h,
                    const uvec &which_event, const uvec &which_right_event,
                    const uvec &which_left, const uvec &which_interval,
                    const bool &any_event, const bool &any_interval,
                    const field<vec> &prior_mean_bs_gammas, field<mat> &prior_Tau_bs_gammas,
                    const vec &tau_bs_gammas,
                    const vec &prior_mean_gammas, mat &prior_Tau_gammas,
                    const vec &lambda_gammas, const double &tau_gammas, const bool &shrink_gammas,
                    const vec &prior_mean_alphas, mat &prior_Tau_alphas,
                    const vec &lambda_alphas, const double &tau_alphas, const bool &shrink_alphas,
                    vec &logLik_surv, double &denominator_surv, const uword &it,
                    /////
                    const mat &Wlong_H, const mat &Wlong_h, const mat &Wlong_H2,
                    vec &scale_alphas, mat &acceptance_alphas, mat &res_alphas,
                    const bool &recurrent,
                    const vec &frailtyH_sigmaF_alphaF, const vec &frailtyh_sigmaF_alphaF,
                    const vec &alphaF, const vec prior_mean_alphaF,
                    mat &prior_Tau_alphaF, const vec &lambda_alphaF,
                    const double &tau_alphaF, const bool &shrink_alphaF) {
  for (uword i = 0; i < alphas.n_rows; ++i) {
    vec proposed_alphas = propose_norm(alphas, scale_alphas, i);
    vec proposed_WlongH_alphas = Wlong_H * proposed_alphas;
    vec proposed_Wlongh_alphas(Wlong_h.n_rows);
    if (any_event) {
      proposed_Wlongh_alphas = Wlong_h * proposed_alphas;
    }
    vec proposed_WlongH2_alphas(Wlong_H2.n_rows);
    if (any_interval) {
      proposed_WlongH2_alphas = Wlong_H2 * proposed_alphas;
    }
    vec logLik_surv_proposed =
      log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
               WH_gammas, Wh_gammas, WH2_gammas,
               proposed_WlongH_alphas, proposed_Wlongh_alphas, proposed_WlongH2_alphas,
               log_Pwk, log_Pwk2, indFast_H, indFast_h,
               which_event, which_right_event, which_left,
               any_interval, which_interval,
               recurrent, frailtyH_sigmaF_alphaF, frailtyh_sigmaF_alphaF);
    double numerator_surv =
      sum(logLik_surv_proposed) +
      logPrior_surv(bs_gammas, gammas, proposed_alphas, prior_mean_bs_gammas,
                    prior_Tau_bs_gammas, tau_bs_gammas,
                    prior_mean_gammas, prior_Tau_gammas, lambda_gammas, tau_gammas, shrink_gammas,
                    prior_mean_alphas, prior_Tau_alphas, lambda_alphas, tau_alphas, shrink_alphas,
                    recurrent, alphaF, prior_mean_alphaF, prior_Tau_alphaF,
                    lambda_alphaF, tau_alphaF, shrink_alphaF);
    double log_ratio = numerator_surv - denominator_surv;
    if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0.0, 1.0)) {
      alphas = proposed_alphas;
      WlongH_alphas = proposed_WlongH_alphas;
      if (any_event) {
        Wlongh_alphas = proposed_Wlongh_alphas;
      }
      if (any_interval) {
        WlongH2_alphas = proposed_WlongH2_alphas;
      }
      logLik_surv = logLik_surv_proposed;
      denominator_surv = numerator_surv;
      acceptance_alphas.at(it, i) = 1;
    }
    if (it > 19) {
      scale_alphas.at(i) =
        robbins_monro(scale_alphas.at(i),
                      acceptance_alphas.at(it, i), it);
    }
    // store results
    res_alphas.at(it, i) = alphas.at(i);
  }
}

void update_alphaF (const vec &bs_gammas, const vec &gammas, const vec &alphas,
                    const vec &W0H_bs_gammas, const vec &W0h_bs_gammas, const vec &W0H2_bs_gammas,
                    const vec &WH_gammas, const vec &Wh_gammas, const vec &WH2_gammas,
                    const vec &WlongH_alphas, const vec &Wlongh_alphas, const vec &WlongH2_alphas,
                    const vec &log_Pwk, const vec &log_Pwk2,
                    const uvec &indFast_H, const uvec &indFast_h,
                    const uvec &which_event, const uvec &which_right_event,
                    const uvec &which_left, const uvec &which_interval,
                    const bool &any_event, const bool &any_interval,
                    const field<vec> &prior_mean_bs_gammas, field<mat> &prior_Tau_bs_gammas,
                    const vec &tau_bs_gammas,
                    const vec &prior_mean_gammas, mat &prior_Tau_gammas,
                    const vec &lambda_gammas, const double &tau_gammas, const bool &shrink_gammas,
                    const vec &prior_mean_alphas, mat &prior_Tau_alphas,
                    const vec &lambda_alphas, const double &tau_alphas, const bool &shrink_alphas,
                    vec &logLik_surv, double &denominator_surv, const uword &it,
                    const mat &Wlong_H, const mat &Wlong_h, const mat &Wlong_H2,
                    //
                    const bool &recurrent,
                    const uvec &which_term_H, const uvec &which_term_h,
                    const vec &frailty_H, const vec &frailty_h,
                    vec &alphaF, vec &alphaF_H, vec &alphaF_h,
                    vec &scale_alphaF, mat &acceptance_alphaF, mat &res_alphaF,
                    const vec prior_mean_alphaF, mat &prior_Tau_alphaF,
                    const vec &lambda_alphaF, const double &tau_alphaF, 
                    const bool &shrink_alphaF,
                    const vec &sigmaF,
                    vec &frailtyH_sigmaF_alphaF, vec &frailtyh_sigmaF_alphaF) {
  vec proposed_alphaF = propose_norm(alphaF, scale_alphaF, 0);
  vec proposed_alphaF_H(WH_gammas.n_rows, fill::ones);
  vec proposed_alphaF_h(Wh_gammas.n_rows, fill::ones);
  proposed_alphaF_H.rows(which_term_H).fill(proposed_alphaF.at(0));
  proposed_alphaF_h.rows(which_term_h).fill(proposed_alphaF.at(0));
  vec proposed_frailtyH_sigmaF_alphaF(WH_gammas.n_rows, fill::zeros);
  vec proposed_frailtyh_sigmaF_alphaF(which_event.n_rows, fill::zeros);
  proposed_frailtyH_sigmaF_alphaF = frailty_H % proposed_alphaF_H * sigmaF;
  proposed_frailtyh_sigmaF_alphaF = frailty_h.rows(which_event) % proposed_alphaF_h.rows(which_event) * sigmaF;
  vec logLik_surv_proposed =
    log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
             WH_gammas, Wh_gammas, WH2_gammas,
             WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
             log_Pwk, log_Pwk2, indFast_H, indFast_h,
             which_event, which_right_event, which_left,
             any_interval, which_interval,
             recurrent, 
             proposed_frailtyH_sigmaF_alphaF, proposed_frailtyh_sigmaF_alphaF);
  double numerator_surv =
    sum(logLik_surv_proposed) +
    logPrior_surv(bs_gammas, gammas, alphas, prior_mean_bs_gammas,
                  prior_Tau_bs_gammas, tau_bs_gammas,
                  prior_mean_gammas, prior_Tau_gammas, lambda_gammas, tau_gammas, shrink_gammas,
                  prior_mean_alphas, prior_Tau_alphas, lambda_alphas, tau_alphas, shrink_alphas,
                  recurrent, proposed_alphaF, prior_mean_alphaF, prior_Tau_alphaF,
                  lambda_alphaF, tau_alphaF, shrink_alphaF);
  double log_ratio = numerator_surv - denominator_surv;
  if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0.0, 1.0)) {
    alphaF = proposed_alphaF;
    alphaF_H = proposed_alphaF_H;
    alphaF_h = proposed_alphaF_h;
    logLik_surv = logLik_surv_proposed;
    denominator_surv = numerator_surv;
    acceptance_alphaF.row(it) = 1;
    frailtyH_sigmaF_alphaF = proposed_frailtyH_sigmaF_alphaF;
    frailtyh_sigmaF_alphaF = proposed_frailtyh_sigmaF_alphaF;
  }
  if (it > 19) {
    scale_alphaF.at(0) =
      robbins_monro(scale_alphaF.at(0),
                    acceptance_alphaF.at(it, 0), it);
    
  }
  // store results
  res_alphaF.at(it, 0) = alphaF.at(0);
}

#endif

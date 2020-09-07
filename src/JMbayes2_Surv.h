#ifndef JMBAYES2SURV_H
#define JMBAYES2SURV_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_bs_gammas (vec &bs_gammas, const vec &gammas, const vec &alphas,
                       vec &W0H_bs_gammas, vec &W0h_bs_gammas, vec &W0H2_bs_gammas,
                       const vec &WH_gammas, const vec &Wh_gammas, const vec &WH2_gammas,
                       const vec &WlongH_alphas, const vec &Wlongh_alphas, const vec &WlongH2_alphas,
                       const vec &log_Pwk, const vec &log_Pwk2, const uvec &id_H,
                       const uvec &which_event, const uvec &which_right_event,
                       const uvec &which_left, const uvec &which_interval,
                       const bool &any_event, const bool &any_interval,
                       const vec &prior_mean_bs_gammas, const mat &prior_Tau_bs_gammas,
                       const double &tau_bs_gammas,
                       const vec &prior_mean_gammas, const mat &prior_Tau_gammas,
                       const vec &prior_mean_alphas, const mat &prior_Tau_alphas,
                       double &denominator_surv, const uword &it,
                       /////
                       const mat &W0_H, const mat &W0_h, const mat &W0_H2,
                       vec &scale_bs_gammas, mat &acceptance_bs_gammas,
                       mat &res_bs_gammas) {
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
    double numerator_surv =
      log_surv(proposed_W0H_bs_gammas, proposed_W0h_bs_gammas, proposed_W0H2_bs_gammas,
               WH_gammas, Wh_gammas, WH2_gammas,
               WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
               log_Pwk, log_Pwk2, id_H,
               which_event, which_right_event, which_left,
               any_interval, which_interval) +
        logPrior(proposed_bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
        logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
        logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    double log_ratio = numerator_surv - denominator_surv;
    if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
      bs_gammas = proposed_bs_gammas;
      W0H_bs_gammas = proposed_W0H_bs_gammas;
      if (any_event) {
        W0h_bs_gammas = proposed_W0h_bs_gammas;
      }
      if (any_interval) {
        W0H2_bs_gammas = proposed_W0H2_bs_gammas;
      }
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
                    const vec &log_Pwk, const vec &log_Pwk2, const uvec &id_H,
                    const uvec &which_event, const uvec &which_right_event,
                    const uvec &which_left, const uvec &which_interval,
                    const bool &any_event, const bool &any_interval,
                    const vec &prior_mean_bs_gammas, const mat &prior_Tau_bs_gammas,
                    const double &tau_bs_gammas,
                    const vec &prior_mean_gammas, const mat &prior_Tau_gammas,
                    const vec &prior_mean_alphas, const mat &prior_Tau_alphas,
                    double &denominator_surv, const uword &it,
                    /////
                    const mat &W_H, const mat &W_h, const mat &W_H2,
                    vec &scale_gammas, mat &acceptance_gammas, mat &res_gammas) {
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
    double numerator_surv =
      log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
               proposed_WH_gammas, proposed_Wh_gammas, proposed_WH2_gammas,
               WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
               log_Pwk, log_Pwk2, id_H,
               which_event, which_right_event, which_left,
               any_interval, which_interval) +
          logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
          logPrior(proposed_gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
          logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    double log_ratio = numerator_surv - denominator_surv;
    if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
      gammas = proposed_gammas;
      WH_gammas = proposed_WH_gammas;
      if (any_event) {
        Wh_gammas = proposed_Wh_gammas;
      }
      if (any_interval) {
        WH2_gammas = proposed_WH2_gammas;
      }
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
                    const vec &log_Pwk, const vec &log_Pwk2, const uvec &id_H,
                    const uvec &which_event, const uvec &which_right_event,
                    const uvec &which_left, const uvec &which_interval,
                    const bool &any_event, const bool &any_interval,
                    const vec &prior_mean_bs_gammas, const mat &prior_Tau_bs_gammas,
                    const double &tau_bs_gammas,
                    const vec &prior_mean_gammas, const mat &prior_Tau_gammas,
                    const vec &prior_mean_alphas, const mat &prior_Tau_alphas,
                    double &denominator_surv, const uword &it,
                    /////
                    const mat &Wlong_H, const mat &Wlong_h, const mat &Wlong_H2,
                    vec &scale_alphas, mat &acceptance_alphas, mat &res_alphas) {
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
    double numerator_surv =
      log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
               WH_gammas, Wh_gammas, WH2_gammas,
               proposed_WlongH_alphas, proposed_Wlongh_alphas, proposed_WlongH2_alphas,
               log_Pwk, log_Pwk2, id_H,
               which_event, which_right_event, which_left,
               any_interval, which_interval) +
        logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
        logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
        logPrior(proposed_alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    double log_ratio = numerator_surv - denominator_surv;
    if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
      alphas = proposed_alphas;
      WlongH_alphas = proposed_WlongH_alphas;
      if (any_event) {
        Wlongh_alphas = proposed_Wlongh_alphas;
      }
      if (any_interval) {
        WlongH2_alphas = proposed_WlongH2_alphas;
      }
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

#endif

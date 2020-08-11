#ifndef JMBAYES2BLOCKS
#define JMBAYES2BLOCKS

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_D (mat &L, vec &sds, const mat &b,
               const uvec &upper_part,
               const double &prior_D_sds_df,
               const double &prior_D_sds_sigma,
               const double &prior_D_L_etaLKJ,
               const int &it, const bool &MALA,
               mat &res_sds, mat &res_L,
               vec &scale_sds, vec &scale_L,
               mat &acceptance_sds, mat &acceptance_L) {
  uword n_sds = sds.n_rows;
  uword n_L = upper_part.n_rows;
  double denominator_sds = logPC_D_sds(sds, L, b, prior_D_sds_df,
                                       prior_D_sds_sigma);
  for (uword i = 0; i < n_sds; ++i) {
    double SS = 0.5 * pow(scale_sds.at(i), 2.0);
    double log_mu_current = log(sds.at(i)) - SS;
    vec proposed_sds = propose_lnorm(sds, log_mu_current, scale_sds, i);
    double numerator_sds = logPC_D_sds(proposed_sds, L, b,
                                       prior_D_sds_df, prior_D_sds_sigma);
    double log_mu_proposed = log(proposed_sds.at(i)) - SS;
    double log_ratio_sds = numerator_sds - denominator_sds +
      R::dlnorm(sds.at(i), log_mu_proposed, scale_sds.at(i), true) -
      R::dlnorm(proposed_sds.at(i), log_mu_current, scale_sds.at(i), true);
    if (std::isfinite(log_ratio_sds) && exp(log_ratio_sds) > R::runif(0.0, 1.0)) {
      sds = proposed_sds;
      denominator_sds = numerator_sds;
      acceptance_sds.at(it, i) = 1;
    }
    if (it > 19) {
      scale_sds.at(i) =
        robbins_monro(scale_sds.at(i), acceptance_sds.at(it, i),
                      it);
    }
    res_sds.at(it, i) = sds.at(i);
  }
  double denominator_L = logPC_D_L(L, sds, b, prior_D_L_etaLKJ);
  for (uword i = 0; i < n_L; ++i) {
    uword upper_part_i = upper_part.at(i);
    double deriv_current(0.0);
    double mu_current(0.0);
    mat proposed_L = L;
    if (MALA) {
      deriv_current = deriv_L(L, sds, b, denominator_L, i, upper_part,
                              prior_D_L_etaLKJ);
      mu_current = L.at(upper_part_i) + 0.5 * scale_L.at(i) * deriv_current;
      proposed_L = propose_L(L, scale_L, upper_part, deriv_current, i, true);
    } else {
      proposed_L = propose_L(L, scale_L, upper_part, deriv_current, i);
    }
    double numerator_L(0.0);
    double deriv_proposed(0.0);
    double mu_proposed(0.0);
    double log_ratio_L(0.0);
    if (proposed_L.is_finite()) {
      numerator_L = logPC_D_L(proposed_L, sds, b, prior_D_L_etaLKJ);
      if (MALA) {
        deriv_proposed = deriv_L(proposed_L, sds, b, numerator_L,
                                 i, upper_part, prior_D_L_etaLKJ);
        mu_proposed = proposed_L.at(upper_part_i) +
          0.5 * scale_L.at(i) * deriv_proposed;
        log_ratio_L = numerator_L - denominator_L +
          log_normpdf(L.at(upper_part_i), mu_proposed, sqrt(scale_L.at(i))) -
          log_normpdf(proposed_L.at(upper_part_i), mu_current, sqrt(scale_L.at(i)));
      } else {
        log_ratio_L = numerator_L - denominator_L;
      }
    }
    if (std::isfinite(log_ratio_L) && exp(log_ratio_L) > R::runif(0.0, 1.0)) {
      L = proposed_L;
      denominator_L = numerator_L;
      acceptance_L.at(it, i) = 1;
    }
    if (it > 19) {
      scale_L.at(i) =
        robbins_monro(scale_L.at(i), acceptance_L.at(it, i),
                      it);
    }
    res_L.at(it, i) = L.at(upper_part_i);
  }
}

void update_bs_gammas (vec& bs_gammas, vec& gammas, vec& alphas,
                       vec& W0H_bs_gammas, vec& W0h_bs_gammas, vec& W0H2_bs_gammas,
                       vec& WH_gammas, vec& Wh_gammas, vec& WH2_gammas,
                       vec& WlongH_alphas, vec& Wlongh_alphas, vec& WlongH2_alphas,
                       vec& log_Pwk, vec& log_Pwk2, uvec& id_H,
                       uvec& which_event, uvec& which_right_event, uvec& which_left, uvec& which_interval,
                       bool& any_event, bool& any_interval,
                       vec& prior_mean_bs_gammas, mat& prior_Tau_bs_gammas, double& tau_bs_gammas,
                       vec& prior_mean_gammas, mat& prior_Tau_gammas,
                       vec& prior_mean_alphas, mat& prior_Tau_alphas,
                       double& denominator_surv, int& it,
                       /////
                       mat& W0_H, mat& W0_h, mat& W0_H2,
                       vec& scale_bs_gammas,
                       mat& acceptance_bs_gammas,
                       mat& res_bs_gammas) {
  for (unsigned int i = 0; i < bs_gammas.n_rows; ++i) {
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
      log_density_surv(proposed_W0H_bs_gammas, proposed_W0h_bs_gammas, proposed_W0H2_bs_gammas,
                       WH_gammas, Wh_gammas, WH2_gammas,
                       WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                       log_Pwk, log_Pwk2, id_H,
                       which_event, which_right_event, which_left,
                       any_interval, which_interval) +
                         logPrior(proposed_bs_gammas, prior_mean_bs_gammas,
                                  prior_Tau_bs_gammas, tau_bs_gammas) +
                                    logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
                                    logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    double log_ratio = numerator_surv - denominator_surv;
    if (is_finite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
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

void update_gammas (vec& bs_gammas, vec& gammas, vec& alphas,
                    vec& W0H_bs_gammas, vec& W0h_bs_gammas, vec& W0H2_bs_gammas,
                    vec& WH_gammas, vec& Wh_gammas, vec& WH2_gammas,
                    vec& WlongH_alphas, vec& Wlongh_alphas, vec& WlongH2_alphas,
                    vec& log_Pwk, vec& log_Pwk2, uvec& id_H,
                    uvec& which_event, uvec& which_right_event, uvec& which_left, uvec& which_interval,
                    bool& any_event, bool& any_interval,
                    vec& prior_mean_bs_gammas, mat& prior_Tau_bs_gammas, double& tau_bs_gammas,
                    vec& prior_mean_gammas, mat& prior_Tau_gammas,
                    vec& prior_mean_alphas, mat& prior_Tau_alphas,
                    double& denominator_surv, int& it,
                    /////
                    mat& W_H, mat& W_h, mat& W_H2,
                    vec& scale_gammas,
                    mat& acceptance_gammas,
                    mat& res_gammas) {
  for (unsigned int i = 0; i < gammas.n_rows; ++i) {
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
      log_density_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                       proposed_WH_gammas, proposed_Wh_gammas, proposed_WH2_gammas,
                       WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                       log_Pwk, log_Pwk2, id_H,
                       which_event, which_right_event, which_left,
                       any_interval, which_interval) +
                         logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas,
                                  tau_bs_gammas) +
                                    logPrior(proposed_gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
                                    logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    double log_ratio = numerator_surv - denominator_surv;
    if (is_finite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
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

void update_alphas (vec& bs_gammas, vec& gammas, vec& alphas,
                    vec& W0H_bs_gammas, vec& W0h_bs_gammas, vec& W0H2_bs_gammas,
                    vec& WH_gammas, vec& Wh_gammas, vec& WH2_gammas,
                    vec& WlongH_alphas, vec& Wlongh_alphas, vec& WlongH2_alphas,
                    vec& log_Pwk, vec& log_Pwk2, uvec& id_H,
                    uvec& which_event, uvec& which_right_event, uvec& which_left, uvec& which_interval,
                    bool& any_event, bool& any_interval,
                    vec& prior_mean_bs_gammas, mat& prior_Tau_bs_gammas, double& tau_bs_gammas,
                    vec& prior_mean_gammas, mat& prior_Tau_gammas,
                    vec& prior_mean_alphas, mat& prior_Tau_alphas,
                    double& denominator_surv, int& it,
                    /////
                    mat& Wlong_H, mat& Wlong_h, mat& Wlong_H2,
                    vec& scale_alphas,
                    mat& acceptance_alphas,
                    mat& res_alphas) {
  for (unsigned int i = 0; i < alphas.n_rows; ++i) {
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
      log_density_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                       WH_gammas, Wh_gammas, WH2_gammas,
                       proposed_WlongH_alphas, proposed_Wlongh_alphas, proposed_WlongH2_alphas,
                       log_Pwk, log_Pwk2, id_H,
                       which_event, which_right_event, which_left,
                       any_interval, which_interval) +
                         logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas,
                                  tau_bs_gammas) +
                                    logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
                                    logPrior(proposed_alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    double log_ratio = numerator_surv - denominator_surv;
    if (is_finite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
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

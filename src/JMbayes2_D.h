#ifndef JMBAYES2D_H
#define JMBAYES2D_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

inline void logPrior_D_sds(const vec &sigmas, const vec &D_sds_sigmas,
                           const double &D_sds_df, const vec &D_sds_mean,
                           const double &D_sds_shape,
                           const bool gamma_prior, vec &out) {
    if (gamma_prior) {
        log_dgamma_void(sigmas, D_sds_shape, D_sds_mean / D_sds_shape, out);
    } else {
        log_dht_void(sigmas, D_sds_sigmas, D_sds_df, out);
    }
}

double logPrior_LKJ (const mat &L, const double &D_L_etaLKJ) {
  uword p = L.n_rows;
  double out(0.0);
  for (uword i = 1; i < p; ++i) {
    out += (p - i - 1.0 + 2.0 * D_L_etaLKJ - 2.0) * std::log(L.at(i, i));
  }
  return out;
}

double logPC_D_L (const mat &L, const vec &sds, const mat &b,
                  const double &D_L_etaLKJ) {
  uword p = L.n_rows;
  mat chol_Sigma = L.each_row() % sds.t();
  double log_p_b = sum(log_dmvnrm_chol(b, chol_Sigma));
  double log_p_L(0.0);
  for (unsigned i = 1; i < p; ++i) {
    log_p_L += (p - i - 1.0 + 2.0 * D_L_etaLKJ - 2.0) * log(L.at(i, i));
  }
  double out = log_p_b + log_p_L;
  return out;
}

double deriv_L (const mat &L, const vec &sds, const mat &b,
                const double &log_target, const uword &i,
                const uvec &upper_part,
                const double &D_L_etaLKJ,
                const char &direction = 'b', const double &eps = 1e-06) {
  uword n = L.n_rows;
  uword upper_part_i = upper_part.at(i);
  uword column = upper_part_i / n;
  mat L_eps = L;
  if (direction == 'f') {
    L_eps(upper_part_i) += L_eps(upper_part_i) * eps;
  } else {
    L_eps(upper_part_i) -= L_eps(upper_part_i) * eps;
  }
  vec ll = L_eps.submat(0, column, column - 1, column);
  double ss = dot(ll, ll);
  if (ss > 1) return datum::nan;
  L_eps.at(column, column) = sqrt(1 - ss);
  double out(0.0);
  if (direction == 'f') {
    out = (logPC_D_L(L_eps, sds, b, D_L_etaLKJ) - log_target) / eps;
  } else {
    out = (log_target - logPC_D_L(L_eps, sds, b, D_L_etaLKJ)) / eps;
  }
  return out;
}

mat propose_L (const mat &L, const vec &scale, const uvec &upper_part,
               const double &deriv, const uword &i, const umat &ind_zero_D,
               const bool &mala = false) {
    mat proposed_L(size(L), fill::zeros);
    vec l = L(upper_part);
    vec proposed_l;
    if (mala) {
        if (std::isfinite(deriv)) {
            proposed_l = propose_norm_mala(l, scale, deriv, i);
        } else {
            proposed_L.fill(datum::nan);
            return proposed_L;
        }
    } else {
        proposed_l = propose_unif(l, scale, i);
    }
    proposed_L(upper_part) = proposed_l;
    uword n = L.n_rows;
    for (uword j = 0; j < n; ++j) {
        auto ll = proposed_L.col(j);
        proposed_L.at(j, j) = std::sqrt(1.0 - arma::dot(ll, ll));
    }
    uword nn = ind_zero_D.n_rows;
    for (uword j = 0; j < nn; ++j) {
        uword j0 = ind_zero_D.at(j, 0);
        uword j1 = ind_zero_D.at(j, 1);
        auto col_j0 = proposed_L.col(j0);
        auto col_j1 = proposed_L.col(j1);
        proposed_L.at(j0, j1) = -arma::dot(col_j0, col_j1) / proposed_L.at(j0, j0);
        auto ll = proposed_L.col(j1).subvec(0, j1 - 1);
        double ss = arma::dot(ll, ll);
        if (ss > 1.0) {
            proposed_L.fill(datum::nan);
            return proposed_L;
        }
        proposed_L.at(j1, j1) = std::sqrt(1.0 - ss);
    }
    return proposed_L;
}

inline void propose_L_void (const mat &L, const vec &scale,
                            const uvec &upper_part, const uword i,
                            const umat &ind_zero_D, mat &proposed_L,
                            vec &proposed_l) {
    proposed_L.zeros();
    proposed_l = L(upper_part);
    double rr = Const_Unif_Proposal * scale.at(i);
    proposed_l.at(i) = R::runif(proposed_l.at(i) - rr, proposed_l.at(i) + rr);
    proposed_L(upper_part) = proposed_l;
    uword n = L.n_rows;
    for (uword j = 0; j < n; ++j) {
        auto ll = proposed_L.col(j);
        proposed_L.at(j, j) = std::sqrt(1.0 - arma::dot(ll, ll));
    }
    uword nn = ind_zero_D.n_rows;
    for (uword j = 0; j < nn; ++j) {
        uword j0 = ind_zero_D.at(j, 0);
        uword j1 = ind_zero_D.at(j, 1);
        auto col_j0 = proposed_L.col(j0);
        auto col_j1 = proposed_L.col(j1);
        proposed_L.at(j0, j1) = -arma::dot(col_j0, col_j1) / proposed_L.at(j0, j0);
        auto ll = proposed_L.col(j1).subvec(0, j1 - 1);
        double ss = arma::dot(ll, ll);
        if (ss > 1.0) {
            proposed_L.fill(datum::nan);
        }
        proposed_L.at(j1, j1) = std::sqrt(1.0 - ss);
    }
}

void update_D (mat &L, vec &sds, const mat &b,
               const uvec &upper_part,
               const double &D_sds_df,
               const vec &D_sds_sigma,
               const double &D_sds_shape,
               const vec &D_sds_mean,
               const bool &gamma_prior,
               const double &D_L_etaLKJ,
               const int &it, const bool &MALA, const umat &ind_zero_D,
               mat &V_R, mat &Z_workspace, mat &Z_transposed,
               mat &B_scaled_workspace, vec &log_prior_sds, vec &proposed_sds,
               mat &proposed_L, vec &proposed_l, vec &logLik_re,
               vec &logLik_re_proposed, mat &res_sds, mat &res_L,
               vec &scale_sds, vec &scale_L, mat &acceptance_sds,
               mat &acceptance_L) {
  uword n_sds = sds.n_rows;
  uword n_L = upper_part.n_rows;
  inv(V_R, trimatu(L));
  double log_det_V_R = -arma::sum(arma::log(L.diag()));
  logPrior_D_sds(sds, D_sds_sigma, D_sds_df, D_sds_mean, D_sds_shape,
                 gamma_prior, log_prior_sds);
  double denominator_sds = arma::accu(logLik_re) + arma::accu(log_prior_sds);
  for (uword i = 0; i < n_sds; ++i) {
    double val = scale_sds.at(i);
    double SS = 0.5 * val * val;
    double current_sd_i = sds.at(i);
    double log_mu_current = std::log(current_sd_i) - SS;
    double proposed_sd_i = R::rlnorm(log_mu_current, scale_sds.at(i));
    sds.at(i) = proposed_sd_i;
    log_re_onlySDS(b, V_R, log_det_V_R, sds, B_scaled_workspace,
                   Z_workspace, logLik_re_proposed);
    logPrior_D_sds(sds, D_sds_sigma, D_sds_df, D_sds_mean,
                   D_sds_shape, gamma_prior, log_prior_sds);
    double numerator_sds = arma::accu(logLik_re_proposed) +
        arma::accu(log_prior_sds);
    double log_mu_proposed = std::log(proposed_sd_i) - SS;
    double log_ratio_sds = numerator_sds - denominator_sds +
        log_dlnorm(current_sd_i, log_mu_proposed, scale_sds.at(i)) -
        log_dlnorm(proposed_sd_i, log_mu_current, scale_sds.at(i));
    if (std::isfinite(log_ratio_sds) &&
        log_ratio_sds > std::log(R::unif_rand())) {
      logLik_re = logLik_re_proposed;
      denominator_sds = numerator_sds;
      acceptance_sds.at(it, i) = 1;
    } else {
        sds.at(i) = current_sd_i;
    }
    if (it > 119) {
      scale_sds.at(i) =
        robbins_monro(scale_sds.at(i), acceptance_sds.at(it, i), it - 100);
    }
    res_sds.at(it, i) = sds.at(i);
  }
  double denominator_L = sum(logLik_re) + logPrior_LKJ(L, D_L_etaLKJ);
  uword n_rows_b = b.n_rows;
  const double* sds_ptr = sds.memptr();
  for (uword c = 0; c < n_sds; ++c) {
      double inv_sd = 1.0 / sds_ptr[c];
      const double* b_col = b.colptr(c);
      double* B_ws_col = B_scaled_workspace.colptr(c);
      for (uword r = 0; r < n_rows_b; ++r) {
          B_ws_col[r] = b_col[r] * inv_sd;
      }
  }
  double sum_log_sds = arma::accu(arma::log(sds));
  for (uword i = 0; i < n_L; ++i) {
    uword upper_part_i = upper_part.at(i);
    double deriv_current(0.0);
    double mu_current(0.0);
    proposed_L = L;
    if (MALA) {
      deriv_current = deriv_L(L, sds, b, denominator_L, i, upper_part,
                              D_L_etaLKJ);
      mu_current = L.at(upper_part_i) + 0.5 * scale_L.at(i) * deriv_current;
      proposed_L = propose_L(L, scale_L, upper_part, deriv_current, i,
                             ind_zero_D, true);
    } else {
      propose_L_void(L, scale_L, upper_part, i, ind_zero_D, proposed_L,
                     proposed_l);
    }
    logLik_re_proposed = logLik_re;
    double numerator_L(0.0);
    double deriv_proposed(0.0);
    double mu_proposed(0.0);
    double log_ratio_L(0.0);
    bool finite_L = proposed_L.is_finite();
    if (finite_L) {
      log_re_onlyL(B_scaled_workspace, proposed_L, sum_log_sds, Z_transposed,
                   logLik_re_proposed);
      numerator_L = sum(logLik_re_proposed) +
        logPrior_LKJ(proposed_L, D_L_etaLKJ);
      if (MALA) {
        deriv_proposed = deriv_L(proposed_L, sds, b, numerator_L,
                                 i, upper_part, D_L_etaLKJ);
        mu_proposed = proposed_L.at(upper_part_i) +
          0.5 * scale_L.at(i) * deriv_proposed;
        log_ratio_L = numerator_L - denominator_L +
          log_normpdf(L.at(upper_part_i), mu_proposed, sqrt(scale_L.at(i))) -
          log_normpdf(proposed_L.at(upper_part_i), mu_current, sqrt(scale_L.at(i)));
      } else {
        log_ratio_L = numerator_L - denominator_L;
      }
    }
    if (finite_L && std::isfinite(log_ratio_L) &&
        log_ratio_L > std::log(R::unif_rand())) {
      L = proposed_L;
      logLik_re = logLik_re_proposed;
      denominator_L = numerator_L;
      acceptance_L.at(it, i) = 1;
    }
    if (it > 119) {
      scale_L.at(i) =
        robbins_monro(scale_L.at(i), acceptance_L.at(it, i), it - 100);
    }
    res_L.at(it, i) = L.at(upper_part_i);
  }
}

#endif

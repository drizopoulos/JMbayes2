#ifndef JMBAYES2D_H
#define JMBAYES2D_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

vec logPrior_D_sds(const vec &sigmas, const vec &D_sds_sigmas, const double &D_sds_df,
                   const vec &D_sds_mean, const double &D_sds_shape,
                   const bool &gamma_prior = true) {
  vec out(sigmas.n_rows);
  if (gamma_prior) {
    out = log_dgamma(sigmas, D_sds_shape, D_sds_mean / D_sds_shape);
  } else {
    out = log_dht(sigmas, D_sds_sigmas, D_sds_df);
  }
  return out;
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
               const double &deriv, const uword &i, const bool &mala = false) {
  mat proposed_L = L;
  vec l = L(upper_part);
  vec proposed_l(l.n_rows);
  if (mala) {
    if (std::isfinite(deriv)) {
      proposed_l = propose_norm_mala(l, scale, deriv, i);
    } else {
      return proposed_L.fill(datum::nan);
    }
  } else {
    proposed_l = propose_unif(l, scale, i);
  }
  proposed_L(upper_part) = proposed_l;
  uword n = L.n_rows;
  uword column = upper_part.at(i) / n;
  vec ll = proposed_L.submat(0, column, column - 1, column);
  double ss = dot(ll, ll);
  if (ss > 1) return proposed_L.fill(datum::nan);
  proposed_L.at(column, column) = sqrt(1 - ss);
  return proposed_L;
}

void update_D (mat &L, vec &sds, const mat &b,
               const uvec &upper_part,
               const double &D_sds_df,
               const vec &D_sds_sigma,
               const double &D_sds_shape,
               const vec &D_sds_mean,
               const bool &gamma_prior,
               const double &D_L_etaLKJ,
               const int &it, const bool &MALA,
               vec &logLik_re,
               mat &res_sds, mat &res_L,
               vec &scale_sds, vec &scale_L,
               mat &acceptance_sds, mat &acceptance_L) {
  uword n_sds = sds.n_rows;
  uword n_L = upper_part.n_rows;
  double denominator_sds = sum(logLik_re) +
    sum(logPrior_D_sds(sds, D_sds_sigma, D_sds_df, D_sds_mean, D_sds_shape, gamma_prior));
  for (uword i = 0; i < n_sds; ++i) {
    double SS = 0.5 * pow(scale_sds.at(i), 2.0);
    double log_mu_current = log(sds.at(i)) - SS;
    vec proposed_sds = propose_lnorm(sds, log_mu_current, scale_sds, i);
    vec logLik_re_proposed = log_re(b, L, proposed_sds);
    double numerator_sds = sum(logLik_re_proposed) +
      sum(logPrior_D_sds(proposed_sds, D_sds_sigma, D_sds_df, D_sds_mean, D_sds_shape, gamma_prior));
    double log_mu_proposed = log(proposed_sds.at(i)) - SS;
    double log_ratio_sds = numerator_sds - denominator_sds +
      R::dlnorm(sds.at(i), log_mu_proposed, scale_sds.at(i), true) -
      R::dlnorm(proposed_sds.at(i), log_mu_current, scale_sds.at(i), true);
    if (std::isfinite(log_ratio_sds) && exp(log_ratio_sds) > R::runif(0.0, 1.0)) {
      sds = proposed_sds;
      logLik_re = logLik_re_proposed;
      denominator_sds = numerator_sds;
      acceptance_sds.at(it, i) = 1;
    }
    if (it > 119) {
      scale_sds.at(i) =
        robbins_monro(scale_sds.at(i), acceptance_sds.at(it, i), it - 100);
    }
    res_sds.at(it, i) = sds.at(i);
  }
  double denominator_L = sum(logLik_re) + logPrior_LKJ(L, D_L_etaLKJ);
  for (uword i = 0; i < n_L; ++i) {
    uword upper_part_i = upper_part.at(i);
    double deriv_current(0.0);
    double mu_current(0.0);
    mat proposed_L = L;
    if (MALA) {
      deriv_current = deriv_L(L, sds, b, denominator_L, i, upper_part,
                              D_L_etaLKJ);
      mu_current = L.at(upper_part_i) + 0.5 * scale_L.at(i) * deriv_current;
      proposed_L = propose_L(L, scale_L, upper_part, deriv_current, i, true);
    } else {
      proposed_L = propose_L(L, scale_L, upper_part, deriv_current, i);
    }
    vec logLik_re_proposed = logLik_re;
    double numerator_L(0.0);
    double deriv_proposed(0.0);
    double mu_proposed(0.0);
    double log_ratio_L(0.0);
    bool finite_L = proposed_L.is_finite();
    if (finite_L) {
      logLik_re_proposed = log_re(b, proposed_L, sds);
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
        exp(log_ratio_L) > R::runif(0.0, 1.0)) {
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

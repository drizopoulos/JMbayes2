#ifndef JMBAYES2SIGMAS_H
#define JMBAYES2SIGMAS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


void update_sigmas (vec &sigmas, const uvec &has_sigmas,
                    const field<mat> &y, const field<vec> &eta,
                    const vec &extra_parms, const CharacterVector &families,
                    const CharacterVector &links, const field<uvec> &idFast,
                    const field<uvec> &unq_ids, const double &prior_sigmas_df,
                    const double &prior_sigmas_sigma, const int &it,
                    vec &logLik_long, mat &res_sigmas, vec &scale_sigmas,
                    mat &acceptance_sigmas) {
  uword n_sigmas = sigmas.n_rows;
  double denominator = sum(logLik_long) +
    sum(log_dht(sigmas, prior_sigmas_sigma, prior_sigmas_df));
  for (uword i = 0; i < n_sigmas; ++i) {
    double SS = 0.5 * std::pow(scale_sigmas.at(i), 2.0);
    double log_mu_current = std::log(sigmas.at(i)) - SS;
    vec proposed_sigmas = propose_lnorm(sigmas, log_mu_current, scale_sigmas, i);
    vec logLik_long_proposed =
      log_long(y, eta, proposed_sigmas, extra_parms, families, links, idFast,
               unq_ids);
    double numerator = sum(logLik_long_proposed) +
      sum(log_dht(proposed_sigmas, prior_sigmas_sigma, prior_sigmas_df));
    double log_mu_proposed = std::log(proposed_sigmas.at(i)) - SS;
    double log_ratio = numerator - denominator +
      R::dlnorm(sigmas.at(i), log_mu_proposed, scale_sigmas.at(i), true) -
      R::dlnorm(proposed_sigmas.at(i), log_mu_current, scale_sigmas.at(i), true);
    if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0.0, 1.0)) {
      sigmas = proposed_sigmas;
      logLik_long = logLik_long_proposed;
      denominator = numerator;
      acceptance_sigmas.at(it, i) = 1;
    }
    if (it > 19) {
      scale_sigmas.at(i) =
        robbins_monro(scale_sigmas.at(i), acceptance_sigmas.at(it, i), it);
    }
    res_sigmas.at(it, i) = sigmas.at(i);
  }
}

#endif

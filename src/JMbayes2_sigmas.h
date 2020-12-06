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
                    const double &sigmas_df, const vec &sigmas_sigma,
                    const uword &it, mat &res_sigmas, vec &scale_sigmas,
                    mat &acceptance_sigmas) {
  uword n_sigmas = sigmas.n_rows;
  for (uword i = 0; i < n_sigmas; ++i) {
    if (!has_sigmas.at(i)) continue;
    vec logLik_long_i =
      log_long_i(y.at(i), eta.at(i), sigmas.at(i), extra_parms.at(i),
                 std::string(families[i]), std::string(links[i]), idFast.at(i));
    double denominator = sum(logLik_long_i) +
      sum(log_dht(sigmas, sigmas_sigma, sigmas_df));
    //
    double SS = 0.5 * std::pow(scale_sigmas.at(i), 2.0);
    double log_mu_current = std::log(sigmas.at(i)) - SS;
    vec proposed_sigmas = propose_lnorm(sigmas, log_mu_current, scale_sigmas, i);
    vec logLik_long_proposed_i =
      log_long_i(y.at(i), eta.at(i), proposed_sigmas.at(i), extra_parms.at(i),
                 std::string(families[i]), std::string(links[i]), idFast.at(i));
    double numerator = sum(logLik_long_proposed_i) +
      sum(log_dht(proposed_sigmas, sigmas_sigma, sigmas_df));
    double log_mu_proposed = std::log(proposed_sigmas.at(i)) - SS;
    double log_ratio = numerator - denominator +
      R::dlnorm(sigmas.at(i), log_mu_proposed, scale_sigmas.at(i), true) -
      R::dlnorm(proposed_sigmas.at(i), log_mu_current, scale_sigmas.at(i), true);
    if (std::isfinite(log_ratio) && std::exp(log_ratio) > R::runif(0.0, 1.0)) {
      sigmas = proposed_sigmas;
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

#ifndef JMBAYES2TAU_H
#define JMBAYES2TAU_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_taus (vec &taus, const vec &thetas, const mat &prior_Tau,
                  const field<uvec> &ind_taus, const uword &it,
                  const vec &post_A_taus, const vec &prior_B_taus,
                  mat &res_taus) {
  uword n_taus = taus.n_rows;
  for (uword j = 0; j < n_taus; ++j) {
    uvec ind_j = ind_taus.at(j);
    vec thetas_j = thetas.rows(ind_j);
    double quad = as_scalar(thetas_j.t() * prior_Tau.submat(ind_j, ind_j) *
                            thetas_j);
    double post_B_tau = prior_B_taus.at(j) + 0.5 * quad;
    taus.at(j) = R::rgamma(post_A_taus.at(j), 1.0 / post_B_tau);
    res_taus.at(it, j) = taus.at(j);
  }
}

#endif

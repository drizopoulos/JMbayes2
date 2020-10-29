#ifndef JMBAYES2PENALTIES_H
#define JMBAYES2PENALTIES_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_penalties (
    const vec &coefs, vec &lambdas, double &tau, vec &nus, double &xi,
    const bool &single, const double &A_lambdas, const double &B_lambdas,
    const double &A_tau, const double &B_tau, const double &A_nus,
    const double &B_nus, const double &A_xi, const double &B_xi) {
  uword p = lambdas.n_rows;
  vec coefs2 = square(coefs);
  if (single) {
    for (uword j = 0; j < p; ++j) {
      double post_A_lambdas = A_lambdas + 0.5;
      double post_B_lambdas = B_lambdas + 0.5 * tau * coefs2.at(j);
      lambdas.at(j) =  R::rgamma(post_A_lambdas, 1.0 / post_B_lambdas);
    }
    double post_A_tau = A_tau + 0.5 * (double)p;
    double post_B_tau = B_tau + 0.5 * sum(lambdas % coefs2);
    tau =  R::rgamma(post_A_tau, 1.0 / post_B_tau);
  } else {
    for (uword j = 0; j < p; ++j) {
      double post_A_lambdas = A_lambdas + 0.5;
      double post_B_lambdas = nus.at(j) + 0.5 * tau * coefs2.at(j);
      lambdas.at(j) =  R::rgamma(post_A_lambdas, 1.0 / post_B_lambdas);
    }
    double post_A_tau = A_tau + 0.5 * (double)p;
    double post_B_tau = xi + 0.5 * sum(lambdas % coefs2);
    tau =  R::rgamma(post_A_tau, 1.0 / post_B_tau);
    //
    for (uword j = 0; j < p; ++j) {
      double post_A_nus = A_nus + 0.5;
      double post_B_nus = B_nus + lambdas.at(j);
      nus.at(j) =  R::rgamma(post_A_nus, 1.0 / post_B_nus);
    }
    double post_A_xi = A_xi + 0.5;
    double post_B_xi = B_xi + tau;
    xi =  R::rgamma(post_A_xi, 1.0 / post_B_xi);
  }
}

#endif

#ifndef JMBAYES2LONG_H
#define JMBAYES2LONG_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_mean_u (field<mat> &mean_u, const field<vec> &betas,
                    const field<mat> &Xbase, const field<uvec> &x_in_z,
                    const field<uvec> &baseline, const field<uvec> &unq_idL) {
  uword n = mean_u.n_elem;
  for (uword i = 0; i < n; ++i) {
    vec betas_i = betas.at(i);
    mat Xbase_i = Xbase.at(i);
    uvec xinz_i = x_in_z.at(i);
    uvec base_i = baseline.at(i);
    uvec rowind_i = unq_idL.at(i);
    uword k = xinz_i.n_rows;
    if (mean_u.at(i).n_cols == k) {
      mean_u.at(i).each_row() = betas_i.rows(xinz_i).t();
    } else {
      mean_u.at(i).cols(0, k - 1).each_row() = betas_i.rows(xinz_i).t();
    }
    if (is_finite(Xbase_i)) {
      mean_u.at(i)(rowind_i) = Xbase_i * betas_i.rows(base_i);
    }
  }
}

#endif

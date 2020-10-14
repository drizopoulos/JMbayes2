#ifndef JMBAYES2LONG_H
#define JMBAYES2LONG_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_Wlong (mat &Wlong_H, mat &Wlong_h, mat &Wlong_H2,
                   const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
                   const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
                   const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
                   const field<vec> &betas, const field<mat> &b,
                   const uvec &id_H, const uvec &id_h,
                   const field<uvec> &FunForms, const field<uvec> &FunForms_ind,
                   const bool &any_event, const bool &any_interval) {
  field<mat> eta_H = linpred_surv(X_H, betas, Z_H, b, id_H);
  Wlong_H = docall_cbindF(create_Wlong(eta_H, FunForms, U_H, FunForms_ind));
  if (any_event) {
    field<mat> eta_h = linpred_surv(X_h, betas, Z_h, b, id_h);
    Wlong_h = docall_cbindF(create_Wlong(eta_h, FunForms, U_h, FunForms_ind));
  }
  if (any_interval) {
    field<mat> eta_H2 = linpred_surv(X_H2, betas, Z_H2, b, id_H);
    Wlong_H2 = docall_cbindF(create_Wlong(eta_H2, FunForms, U_H2, FunForms_ind));
  }
}

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

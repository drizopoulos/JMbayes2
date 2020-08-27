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

#endif

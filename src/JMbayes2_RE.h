#ifndef JMBAYES2SURV_H
#define JMBAYES2SURV_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_b (cube &b, mat &b_mat, const field<mat> &Xbetas, const field<mat> &Z, const field<uvec> &id, 
               const field<mat> &y, const vec &scales, const vec &extra_parms, 
               const CharacterVector &families, const CharacterVector &links, const field<uvec> &ids,
               const field<uvec> &unq_ids, const mat &L) {
  
}
#ifndef JMBAYES2RE_H
#define JMBAYES2RE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_b (field<mat> &b, mat &b_mat, field<vec> &eta,
               vec &logLik_long, vec &logLik_surv, vec &logLik_re,
               mat &Wlong_H, mat &Wlong_h, mat &Wlong_H2,
               vec &WlongH_alphas, vec &Wlongh_alphas, vec &WlongH2_alphas,
               mat &scale_b,
               const field<uvec> &ind_RE,
               const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
               const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
               const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
               const mat &Wlong_bar, const mat &Wlong_sds, const field<vec> &betas,
               const vec &alphas, const uvec &id_H_, const uvec &id_h,
               const field<uvec> &FunForms, const List Funs_FunForms,
               const field<mat> &X, const field<mat> &Z,
               const field<uvec> &idL, const field<mat> &y, const vec &sigmas,
               const vec &extra_parms, const CharacterVector &families,
               const CharacterVector &links, const field<uvec> &ids,
               const field<uvec> &unq_ids, const vec &W0H_bs_gammas,
               const vec &W0h_bs_gammas,
               const vec &W0H2_bs_gammas, const vec &WH_gammas,
               const vec &Wh_gammas, const vec &WH2_gammas,
               const vec &log_Pwk, const vec &log_Pwk2, const vec &log_weights,
               const uvec &id_h2, const uvec &intgr_ind, const bool &intgr,
               const uvec &indFast_H, const uvec &indFast_h, const uvec &which_event,
               const uvec &which_right_event, const uvec &which_left,
               const uvec &which_interval, const bool &any_event,
               const bool &any_interval, const umat &ni_event,
               const mat &L, const vec &sds,
               const uword &it, mat &acceptance_b,
               const uword &n_burnin,
               const uword &GK_k,
               const bool &recurrent,
               const vec &frailtyH_sigmaF_alphaF,
               const vec &frailtyh_sigmaF_alphaF) {
    uword n = b_mat.n_rows;
    uword nRE = b_mat.n_cols;
    mat V_R = inv(trimatu(L));
    mat V_Sigma = V_R.each_col() / sds;
    double log_det_V_Sigma = -arma::sum(arma::log(L.diag())) -
        arma::sum(arma::log(sds));
    double other_terms = -(double)nRE / 2.0 * log2pi + log_det_V_Sigma;
    vec denominator_b = logLik_long + logLik_surv + logLik_re;
    uvec map_o(nRE);
    uvec map_k(nRE);
    for (uword o = 0; o < ind_RE.n_elem; ++o) {
        for (uword k = 0; k < ind_RE.at(o).n_elem; ++k) {
            uword j_col = ind_RE.at(o).at(k);
            map_o[j_col] = o;
            map_k[j_col] = k;
        }
    }
    field<mat> proposed_b(ind_RE.n_elem);
    for (uword i = 0; i < ind_RE.n_elem; ++i) {
        proposed_b.at(i).set_size(b_mat.n_rows, ind_RE.at(i).n_elem);
    }
    for (uword j = 0; j < nRE; ++j) {
        vec old_b_j = b_mat.col(j);
        vec new_b_j = old_b_j + scale_b.col(j) % arma::randn<arma::vec>(n);
        vec delta_b = new_b_j - old_b_j;
        b_mat.col(j) = new_b_j;
        mat2field_inplace(proposed_b, b_mat, ind_RE);
        uword o = map_o[j];
        uword k = map_k[j];
        uword N_o = Z.at(o).n_rows;
        double* eta_ptr = eta.at(o).memptr();
        const double* Z_col = Z.at(o).colptr(k);
        const double* delta_ptr = delta_b.memptr();
        const uword* id_ptr = idL.at(o).memptr();
        for (uword obs = 0; obs < N_o; ++obs) {
            eta_ptr[obs] += Z_col[obs] * delta_ptr[id_ptr[obs]];
        }
        vec logLik_long_proposed =
            log_long(y, eta, sigmas, extra_parms, families, links, ids, unq_ids,
                     n);
        mat Wlong_H_proposed =
            calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, Wlong_sds, betas,
                            proposed_b, id_H_, FunForms, Funs_FunForms);
        vec WlongH_alphas_proposed = Wlong_H_proposed * alphas;
        mat Wlong_h_proposed;
        vec Wlongh_alphas_proposed;
        if (any_event) {
            Wlong_h_proposed =
                calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, Wlong_sds, betas,
                                proposed_b, id_h, FunForms, Funs_FunForms);
            Wlongh_alphas_proposed = Wlong_h_proposed * alphas;
        }
        mat Wlong_H2_proposed;
        vec WlongH2_alphas_proposed;
        if (any_interval) {
            Wlong_H2_proposed =
                calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, Wlong_sds, betas,
                                proposed_b, id_H_, FunForms, Funs_FunForms);
            WlongH2_alphas_proposed = Wlong_H2_proposed * alphas;
        }
        vec logLik_surv_proposed =
            log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                     WH_gammas, Wh_gammas, WH2_gammas,
                     WlongH_alphas_proposed, Wlongh_alphas_proposed,
                     WlongH2_alphas_proposed,
                     log_Pwk, log_Pwk2, log_weights, id_h2, intgr_ind, intgr,
                     indFast_H, indFast_h,
                     which_event, which_right_event, which_left,
                     any_interval, which_interval,
                     recurrent, frailtyH_sigmaF_alphaF, frailtyh_sigmaF_alphaF);
        vec logLik_re_proposed = log_re_onlyRE(b_mat, V_Sigma, other_terms);
        vec numerator_b = logLik_long_proposed + logLik_surv_proposed +
            logLik_re_proposed;
        vec log_ratio = numerator_b - denominator_b;
        uvec accepted(n, arma::fill::zeros);
        uword* acc_ptr = accepted.memptr();
        for (uword i = 0; i < n; ++i) {
            double acc_i = 0.0;
            if (std::isfinite(log_ratio.at(i)) && std::log(R::unif_rand()) < log_ratio.at(i)) {
                acc_i = 1.0;
                acc_ptr[i] = 1;
                if (it > n_burnin - 1) acceptance_b.at(i, j) += 1.0;
                denominator_b.at(i) = numerator_b.at(i);
                logLik_long.at(i) = logLik_long_proposed.at(i);
                logLik_surv.at(i) = logLik_surv_proposed.at(i);
                logLik_re.at(i) = logLik_re_proposed.at(i);
                uword first_H = GK_k * ni_event.at(i, 0);
                uword last_H = GK_k * ni_event.at(i, 1) - 1;
                Wlong_H.rows(first_H, last_H) =
                    Wlong_H_proposed.rows(first_H, last_H);
                WlongH_alphas.rows(first_H, last_H) =
                    WlongH_alphas_proposed.rows(first_H, last_H);
                if (any_event) {
                    uword first_h = ni_event.at(i, 0);
                    uword last_h = ni_event.at(i, 1) - 1;
                    Wlong_h.rows(first_h, last_h) =
                        Wlong_h_proposed.rows(first_h, last_h);
                    Wlongh_alphas.rows(first_h, last_h) =
                        Wlongh_alphas_proposed.rows(first_h, last_h);
                }
                if (any_interval) {
                    Wlong_H2.rows(first_H, last_H) =
                        Wlong_H2_proposed.rows(first_H, last_H);
                    WlongH2_alphas.rows(first_H, last_H) =
                        WlongH2_alphas_proposed.rows(first_H, last_H);
                }
            } else {
                b_mat.at(i, j) = old_b_j[i];
            }
            if (it > 119) {
                scale_b.at(i, j) =
                    robbins_monro(scale_b.at(i, j), acc_i, it - 100);
            }
        }
        for (uword obs = 0; obs < N_o; ++obs) {
            uword subj = id_ptr[obs];
            if (!acc_ptr[subj]) {
                eta_ptr[obs] -= Z_col[obs] * delta_ptr[subj];
            }
        }
    }
    mat2field_inplace(b, b_mat, ind_RE);
}

void update_frailty (vec &frailty, mat &res_frailty, mat &acceptance_frailty,
                     vec &scale_frailty, vec &frailty_H, vec &frailty_h,
                     vec &logLik_surv, vec &logLik_frailty,
                     const bool &recurrent,
                     const vec &alphaF_H, const vec &alphaF_h,
                     const vec &WlongH_alphas, const vec &Wlongh_alphas,
                     const vec &WlongH2_alphas,
                     const uvec &id_H_, const uvec &id_h,
                     const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
                     const vec &W0H2_bs_gammas, const vec &WH_gammas,
                     const vec &Wh_gammas, const vec &WH2_gammas,
                     const vec &log_Pwk, const vec &log_Pwk2, const vec &log_weights,
                     const uvec &id_h2, const uvec &intgr_ind, const bool &intgr,
                     const uvec &indFast_H, const uvec &indFast_h,
                     const uvec &which_event, const uvec &which_right_event,
                     const uvec &which_left, const uvec &which_interval,
                     const bool &any_interval,
                     const uword &n_burnin, const uword &it,
                     const vec &sigmaF,
                     vec &frailtyH_sigmaF_alphaF, vec &frailtyh_sigmaF_alphaF) {
  uword n = frailty.n_rows;
  // calculate denominator
  vec denominator_frailty = logLik_surv + logLik_frailty;
  // propose new frailty
  vec frailty_proposed = propose_rnorm_vec(frailty, scale_frailty);
  // calculate logLik_Surv_proposed
  vec frailty_H_proposed(frailty_H.n_rows, fill::zeros);
  vec frailty_h_proposed(frailty_h.n_rows, fill::zeros);
  frailty_H_proposed = frailty_proposed.rows(id_H_);
  frailty_h_proposed = frailty_proposed.rows(id_h);
  vec proposed_frailtyH_sigmaF_alphaF(WH_gammas.n_rows, fill::zeros);
  vec proposed_frailtyh_sigmaF_alphaF(which_event.n_rows, fill::zeros);
  proposed_frailtyH_sigmaF_alphaF = frailty_H_proposed % alphaF_H * sigmaF;
  proposed_frailtyh_sigmaF_alphaF = frailty_h_proposed.rows(which_event) % alphaF_h.rows(which_event) * sigmaF;
  vec logLik_surv_proposed =
    log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
             WH_gammas, Wh_gammas, WH2_gammas,
             WlongH_alphas, Wlongh_alphas,
             WlongH2_alphas,
             log_Pwk, log_Pwk2, log_weights, id_h2, intgr_ind, intgr,
             indFast_H, indFast_h,
             which_event, which_right_event, which_left,
             any_interval, which_interval,
             recurrent,
             proposed_frailtyH_sigmaF_alphaF, proposed_frailtyh_sigmaF_alphaF);
  // logLik_frailty_proposed
  vec logLik_frailty_proposed = log_dnorm(frailty_proposed, vec(frailty.n_elem, fill::zeros), 1.0);
  // calculate the numerator
  vec numerator_frailty = logLik_surv_proposed + logLik_frailty_proposed;
  // log_ratio
  vec log_ratio = numerator_frailty - denominator_frailty;
  for (uword i = 0; i < n; ++i) {
    double acc_i(0.0);
    if (std::isfinite(log_ratio.at(i)) &&
        exp(log_ratio.at(i)) > R::runif(0.0, 1.0)) {
      acc_i = 1.0;
      if (it > n_burnin - 1) acceptance_frailty.at(i, 0) += 1.0;
      frailty.at(i) = frailty_proposed.at(i);
      logLik_surv.at(i) = logLik_surv_proposed.at(i);
      logLik_frailty.at(i) = logLik_frailty_proposed.at(i);
    }
    if (it > 19) {
      scale_frailty.at(i) =
        robbins_monro(scale_frailty.at(i), acc_i, it);
    }
  }
  frailty_H = frailty.rows(id_H_);
  frailty_h = frailty.rows(id_h);
  res_frailty.row(it) = frailty.t();
  frailtyH_sigmaF_alphaF = frailty_H % alphaF_H * sigmaF;
  frailtyh_sigmaF_alphaF = frailty_h.rows(which_event) % alphaF_h.rows(which_event) * sigmaF;
}

#endif

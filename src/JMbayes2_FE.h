#ifndef JMBAYES2FE_H
#define JMBAYES2FE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
#include "JMbayes2_Funs.h"
#include "JMbayes2_LogDens.h"
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_betas (field<vec> &betas, mat &res_betas, field<vec> &acceptance_betas,
                   field<vec> &scale_betas, field<vec> &eta, vec &logLik_long,
                   field<mat> &eta_H, field<mat> &eta_h, field<mat> &eta_H2,
                   vec &logLik_surv, mat &Wlong_H, mat &Wlong_h, mat &Wlong_H2,
                   vec &WlongH_alphas, vec &Wlongh_alphas, vec &WlongH2_alphas,
                   const vec &Tau_mean_betas_HC, const mat &prior_Tau_betas_HC,
                   mat &b_mat,
                   const mat &L, const vec &sds, const mat &X_dot,
                   const field<uvec> &ind_FE,
                   const field<uvec> &ind_RE,
                   const uvec &ind_FE_HC,
                   const uvec &id_patt,
                   const field<uvec> &ind_RE_patt,
                   const field<uvec> &ind_FE_patt,
                   const uword &it,
                   const uvec &has_tilde_betas,
                   const field<mat> &X,
                   const field<mat> &Z,
                   field<mat> &b,
                   const field<uvec> &idL,
                   const field<mat> &y,
                   const vec &sigmas,
                   const vec &extra_parms,
                   const CharacterVector &families,
                   const CharacterVector &links,
                   const field<uvec> &idL_lp_fast,
                   const field<vec> &prior_mean_betas_nHC,
                   field<mat> &prior_Tau_betas_nHC,
                   const field<uvec> &x_notin_z,
                   const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
                   const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
                   const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
                   const mat &Wlong_bar, const mat &Wlong_sds,
                   const uvec &id_H_, const uvec &id_h,
                   const field<uvec> &FunForms,
                   const List Funs_FunForms,
                   const vec &alphas,
                   const bool &any_event, const bool &any_interval,
                   const vec &W0H_bs_gammas, const vec &W0h_bs_gammas, const vec &W0H2_bs_gammas,
                   const vec &WH_gammas, const vec &Wh_gammas, const vec &WH2_gammas,
                   const vec &log_Pwk, const vec &log_Pwk2, const vec &log_weights,
                   const uvec &id_h2, const uvec &intgr_ind, const bool &intgr,
                   const uvec &id_H_fast, const uvec &id_h_fast,
                   const uvec &which_event, const uvec &which_right_event, const uvec &which_left,
                   const uvec &which_interval, const field<uvec> &unq_idL,
                   const uword &n_burnin,
                   const bool &recurrent,
                   const vec &frailtyH_sigmaF_alphaF, const vec &frailtyh_sigmaF_alphaF,
                   const bool &save_random_effects, cube &res_b, cube &res_b_last,
                   mat &cumsum_b, cube &outprod_b, const uword &n_iter,
                   vec &lambda_H_workspace, vec &H_workspace,
                   vec &lambda_H2_workspace, vec &H2_workspace,
                   vec &surv_out_workspace, vec &logLik_surv_proposed,
                   const field<mat> &X_dots, mat &sum_JXDXJ, vec &sum_JXDu,
                   mat &u_mat, mat &mean_u_mat, mat &mean_u_mat2) {

    uword n_b = b_mat.n_rows;
    uword q = b_mat.n_cols;

    // FE in HC - Gibbs sampling
    vec betas_vec = docall_rbindF(betas);
    vec mean_u = X_dot * betas_vec.rows(ind_FE_HC);

    // 1. VECTORIZED RESHAPING: Replaces the manual u_mat loop
    mean_u_mat = arma::reshape(mean_u, q, n_b).t();
    u_mat = b_mat + mean_u_mat;

    uword patt_count = ind_RE_patt.n_elem;
    uword p_HC = ind_FE_HC.n_elem;

    sum_JXDXJ.zeros();
    sum_JXDu.zeros();
    mat U = L.each_row() % sds.t();

    // 2. ISOLATED D_INV PRE-CALCULATION
    // Store just the upper-triangular Cholesky factor, not the inverse
    field<mat> U_patt_field(patt_count);
    for (uword p = 0; p < patt_count; ++p) {
        if (!ind_RE_patt.at(p).is_empty()) {
            // No inv(), no matrix multiplication!
            U_patt_field.at(p) = trimatu(chol_update(U, ind_RE_patt.at(p)));
        }
    }

    mat X_tilde;
    vec u_tilde;
    uvec absolute_rows;
    vec u_i;

    for (uword i = 0; i < n_b; ++i) {
        uword patt_i = id_patt.at(i);
        if (ind_FE_patt.at(patt_i).is_empty()) continue;

        const uvec& ind_FE_i = ind_FE_patt.at(patt_i);
        const uvec& ind_RE_i = ind_RE_patt.at(patt_i);
        absolute_rows = i * q + ind_RE_i;

        u_i = u_mat.row(i).t();
        u_i = u_i.rows(ind_RE_i);

        // Fetch the upper-triangular Cholesky factor
        const mat& R = U_patt_field.at(patt_i);

        // Fast Triangular Solves (Calculates R^-T * X  and  R^-T * u)
        arma::solve(X_tilde, arma::trimatl(R.t()), X_dots.at(i));
        arma::solve(u_tilde, arma::trimatl(R.t()), u_i);

        // Simple Cross-products
        sum_JXDu.rows(ind_FE_i) += X_tilde.t() * u_tilde;
        sum_JXDXJ.submat(ind_FE_i, ind_FE_i) += X_tilde.t() * X_tilde;
    }

    //mat Sigma_1 = inv_sympd(prior_Tau_betas_HC + sum_JXDXJ);
    //vec mean_1 = Sigma_1 * (Tau_mean_betas_HC + sum_JXDu);
    //betas_vec.rows(ind_FE_HC) = mean_1 + chol(Sigma_1, "lower") * randn<vec>(p_HC);
    // 1. Define the Precision Matrix (Q) and Canonical Mean (b)
    mat Q = prior_Tau_betas_HC + sum_JXDXJ;
    vec b_vec = Tau_mean_betas_HC + sum_JXDu;
    // 2. Calculate the Cholesky factor of the PRECISION matrix EXACTLY ONCE
    mat L_prec = chol(Q, "lower");
    // 3. Forward substitution: Solve L_prec * yy = b_vec  (This is yy = L^-1 * b)
    vec yy = arma::solve(arma::trimatl(L_prec), b_vec);
    // 4. Add the standard normal noise
    vec ww = yy + arma::randn<vec>(p_HC);
    // 5. Back substitution: Solve L_prec^T * x = ww (This is x = L^-T * ww)
    betas_vec.rows(ind_FE_HC) = arma::solve(arma::trimatu(L_prec.t()), ww);
    vec2field_inplace(betas, betas_vec, ind_FE);

    // 3. VECTORIZED b_mat RECALCULATION (Deletes an entire n_b loop)
    mean_u = X_dot * betas_vec.rows(ind_FE_HC);
    mean_u_mat2 = arma::reshape(mean_u, q, n_b).t();
    b_mat = u_mat - mean_u_mat2;

    if (save_random_effects) {
        res_b.slice(it) = b_mat;
    } else if (it > n_burnin - 1) {
        cumsum_b += b_mat;
        for (uword j = 0; j < n_b; j++) {
            outprod_b.slice(j) += b_mat.row(j).t() * b_mat.row(j);
        }
    }

    if (it == n_iter - 1) {
        res_b_last.slice(0) = b_mat;
    }
    mat2field_inplace(b, b_mat, ind_RE);

    // update eta and logLik_surv baselines
    linpred_mixed_inplace(eta, X, betas, Z, b, idL);

    calculate_Wlong_inplace(Wlong_H, eta_H, X_H, Z_H, U_H, Wlong_bar, Wlong_sds,
                            betas, b, id_H_, FunForms, Funs_FunForms);
    WlongH_alphas = Wlong_H * alphas;
    if (any_event) {
        calculate_Wlong_inplace(Wlong_h, eta_h, X_h, Z_h, U_h, Wlong_bar,
                                Wlong_sds, betas, b, id_h, FunForms,
                                Funs_FunForms);
        Wlongh_alphas = Wlong_h * alphas;
    }
    if (any_interval) {
        calculate_Wlong_inplace(Wlong_H2, eta_H2, X_H2, Z_H2, U_H2, Wlong_bar,
                                Wlong_sds, betas, b, id_H_, FunForms,
                                Funs_FunForms);
        WlongH2_alphas = Wlong_H2 * alphas;
    }

    log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas, WH_gammas, Wh_gammas, WH2_gammas,
             WlongH_alphas, Wlongh_alphas, WlongH2_alphas, log_Pwk, log_Pwk2, log_weights,
             id_h2, intgr_ind, intgr, id_H_fast, id_h_fast, which_event, which_right_event,
             which_left, any_interval, which_interval, recurrent, frailtyH_sigmaF_alphaF,
             frailtyh_sigmaF_alphaF, lambda_H_workspace, H_workspace,
             lambda_H2_workspace, H2_workspace, surv_out_workspace, logLik_surv);

    // /////////////////////////////////////////////////////////////////////////////
    // FE outside HC - Metropolis-Hastings sampling
    if (any(has_tilde_betas)) {
        uword n_outcomes = betas.n_elem;
        for (uword j = 0; j < n_outcomes; ++j) {
            if (!has_tilde_betas.at(j)) continue;

            uvec ind_j = x_notin_z.at(j);
            uword n_betas = ind_j.n_rows;

            double sum_logLik_long_j =
                sum(log_long_i(y.at(j), eta.at(j), sigmas.at(j),
                               extra_parms.at(j), std::string(families[j]),
                               std::string(links[j]), idL_lp_fast.at(j)));
            vec ll(n_betas);
            double logPrior_j =
                logPrior(betas.at(j).rows(ind_j), prior_mean_betas_nHC.at(j),
                         prior_Tau_betas_nHC.at(j), ll, 1.0, false);

            double denominator_j = sum_logLik_long_j + sum(logLik_surv) + logPrior_j;

            for (uword i = 0; i < n_betas; ++i) {
                // 4. IN-PLACE SCALAR MUTATION (Destroys field<vec> deep copy bottleneck)
                uword idx = ind_j.at(i);
                double old_beta = betas.at(j).at(idx);
                double diff = scale_betas.at(j).at(i) * R::norm_rand();
                betas.at(j).at(idx) += diff;

                double logPrior_j_prop = logPrior(betas.at(j).rows(ind_j), prior_mean_betas_nHC.at(j),
                                                  prior_Tau_betas_nHC.at(j), ll, 1.0, false);

                // 5. RANK-1 ETA UPDATE (Deletes the massive linpred_mixed_i matrix multiplication)
                vec eta_j_prop = eta.at(j) + X.at(j).col(idx) * diff;

                double sum_logLik_long_j_prop =
                    sum(log_long_i(y.at(j), eta_j_prop, sigmas.at(j),
                                   extra_parms.at(j), std::string(families[j]),
                                   std::string(links[j]), idL_lp_fast.at(j)));

                // 6. DEFERRED MATRIX ALLOCATIONS
                mat Wlong_H_prop =
                    calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, Wlong_sds, betas,
                                    b, id_H_, FunForms, Funs_FunForms);
                vec WlongH_alphas_prop = Wlong_H_prop * alphas;

                mat Wlong_h_prop; vec Wlongh_alphas_prop;
                if (any_event) {
                    Wlong_h_prop =
                        calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, Wlong_sds,
                                        betas, b, id_h, FunForms, Funs_FunForms);
                    Wlongh_alphas_prop = Wlong_h_prop * alphas;
                }

                mat Wlong_H2_prop; vec WlongH2_alphas_prop;
                if (any_interval) {
                    Wlong_H2_prop =
                        calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, Wlong_sds,
                                        betas, b, id_H_, FunForms, Funs_FunForms);
                    WlongH2_alphas_prop = Wlong_H2_prop * alphas;
                }

                log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas, WH_gammas, Wh_gammas, WH2_gammas,
                         WlongH_alphas_prop, Wlongh_alphas_prop, WlongH2_alphas_prop,
                         log_Pwk, log_Pwk2, log_weights, id_h2, intgr_ind, intgr, id_H_fast, id_h_fast,
                         which_event, which_right_event, which_left, any_interval, which_interval,
                         recurrent, frailtyH_sigmaF_alphaF, frailtyh_sigmaF_alphaF,
                         lambda_H_workspace, H_workspace,
                         lambda_H2_workspace, H2_workspace, surv_out_workspace,
                         logLik_surv_proposed);

                double numerator_j =
                    sum_logLik_long_j_prop + sum(logLik_surv_proposed) + logPrior_j_prop;
                double log_ratio_j = numerator_j - denominator_j;
                double acc_i = 0.0;

                // 7. FAST LOG-SPACE ACCEPTANCE
                if (std::isfinite(log_ratio_j) && std::log(R::unif_rand()) < log_ratio_j) {
                    acc_i = 1.0;
                    if (it > n_burnin - 1) acceptance_betas.at(j).at(i) += 1.0;

                    // betas is already mutated
                    eta.at(j) = eta_j_prop;

                    Wlong_H = Wlong_H_prop;
                    WlongH_alphas = WlongH_alphas_prop;
                    if (any_event) {
                        Wlong_h = Wlong_h_prop;
                        Wlongh_alphas = Wlongh_alphas_prop;
                    }
                    if (any_interval) {
                        Wlong_H2 = Wlong_H2_prop;
                        WlongH2_alphas = WlongH2_alphas_prop;
                    }
                    logLik_surv = logLik_surv_proposed;
                    denominator_j = numerator_j;
                } else {
                    // Reject: Simply revert the scalar change in the field
                    betas.at(j).at(idx) = old_beta;
                }

                if (it > 119) {
                    scale_betas.at(j).at(i) =
                        robbins_monro(scale_betas.at(j).at(i), acc_i, it - 100);
                }
            }
        }
    }
    log_long(y, eta, sigmas, extra_parms, families, links, idL_lp_fast,
             unq_idL, logLik_long);
    res_betas.row(it) = docall_rbindF(betas).t();
}

#endif

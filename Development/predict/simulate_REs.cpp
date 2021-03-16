#include <RcppArmadillo.h>
#include "JMbayes2_LogDens.h"
#include "JMbayes2_Funs.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat simulate_REs (List Data, List MCMC, List control) {
    //////////////////////////////
    // Event Process Components //
    //////////////////////////////
    uvec which_event = as<uvec>(Data["which_event"]) - 1;
    uvec which_right = as<uvec>(Data["which_right"]) - 1;
    uvec which_right_event = join_cols(which_event, which_right);
    uvec which_left = as<uvec>(Data["which_left"]) - 1;
    uvec which_interval = as<uvec>(Data["which_interval"]) - 1;
    bool any_event = which_event.n_rows > 0;
    bool any_interval = which_interval.n_rows > 0;
    umat ni_event = as<umat>(Data["ni_event"]);
    uvec idT = as<uvec>(Data["idT"]) - 1;
    vec log_Pwk = as<vec>(Data["log_Pwk"]);
    vec log_Pwk2 = as<vec>(Data["log_Pwk2"]);
    uvec id_H = as<uvec>(Data["id_H"]) - 1;
    uvec id_H_ = as<uvec>(Data["id_H_"]) - 1;
    uvec id_h = as<uvec>(Data["id_h"]) - 1;
    uvec indFast_H = create_fast_ind(id_H + 1);
    uvec indFast_h = create_fast_ind(id_h + 1);
    uword GK_k = as<uword>(control["GK_k"]);
    bool any_gammas = as<bool>(Data["any_gammas"]);
    field<uvec> FunForms = List2Field_uvec(as<List>(Data["FunForms_cpp"]), true);
    field<uvec> FunForms_ind = List2Field_uvec(as<List>(Data["FunForms_ind"]), true);
    List Funs_FunForms = as<List>(Data["Funs_FunForms"]);
    //
    field<uvec> ind_RE = List2Field_uvec(as<List>(Data["ind_RE"]), true);
    mat W0_H = as<mat>(Data["W0_H"]);
    mat W0_h = as<mat>(Data["W0_h"]);
    mat W0_H2 = as<mat>(Data["W0_H2"]);
    mat W_H = as<mat>(Data["W_H"]);
    mat W_h = as<mat>(Data["W_h"]);
    mat W_H2 = as<mat>(Data["W_H2"]);
    field<mat> X_H = List2Field_mat(as<List>(Data["X_H"]));
    field<mat> X_h = List2Field_mat(as<List>(Data["X_h"]));
    field<mat> X_H2 = List2Field_mat(as<List>(Data["X_H2"]));
    field<mat> Z_H = List2Field_mat(as<List>(Data["Z_H"]));
    field<mat> Z_h = List2Field_mat(as<List>(Data["Z_h"]));
    field<mat> Z_H2 = List2Field_mat(as<List>(Data["Z_H2"]));
    field<mat> U_H = List2Field_mat(as<List>(Data["U_H"]));
    field<mat> U_h = List2Field_mat(as<List>(Data["U_h"]));
    field<mat> U_H2 = List2Field_mat(as<List>(Data["U_H2"]));
    mat Wlong_bar = docall_cbindL(as<List>(Data["Wlong_bar"]));
    mat Wlong_sds = docall_cbindL(as<List>(Data["Wlong_sds"]));

    /////////////////////////////////////
    // Longitudinal Process Components //
    /////////////////////////////////////
    field<mat> y = List2Field_mat(as<List>(Data["y"]));
    field<mat> X = List2Field_mat(as<List>(Data["X"]));
    field<mat> Z = List2Field_mat(as<List>(Data["Z"]));
    //
    vec extra_parms = as<vec>(Data["extra_parms"]);
    CharacterVector families = as<CharacterVector>(Data["family_names"]);
    CharacterVector links = as<CharacterVector>(Data["links"]);
    field<uvec> idL = List2Field_uvec(as<List>(Data["idL"]), true);
    field<uvec> unq_idL = List2Field_uvec(as<List>(Data["unq_idL"]), true);
    field<uvec> idL_lp = List2Field_uvec(as<List>(Data["idL_lp"]), true);
    field<uvec> ids(idL_lp.n_elem);
    for (uword i = 0; i < idL_lp.n_elem; ++i) {
        ids.at(i) = create_fast_ind(idL_lp.at(i) + 1);
    }

    ////////////////////////////
    // MCMC Sample Parameters //
    ///////////////////////////
    field<mat> b = List2Field_mat(as<List>(MCMC["b"]));
    mat b_mat = docall_cbindF(b);
    mat bs_gammas = trans(as<mat>(MCMC["bs_gammas"]));
    mat gammas = trans(as<mat>(MCMC["gammas"]));
    mat alphas = trans(as<mat>(MCMC["alphas"]));
    field<mat> betas = List2Field_mat(as<List>(MCMC["betas"]));
    for (uword j = 0; j < betas.n_elem; ++j) betas.at(j) = trans(betas.at(j));
    mat sigmas = trans(as<mat>(MCMC["sigmas"]));
    cube D = as<cube>(MCMC["D"]);
    uword K = D.n_slices;
    uword q = D.slice(0).n_rows;
    mat sds(q, K);
    cube L(q, q, K);
    for (uword i = 0; i < K; ++i) {
        mat D_i = D.slice(i);
        sds.col(i) = sqrt(D_i.diag());
        mat R = cov2cor(D_i);
        L.slice(i) = chol(R);
    }
    //////////////////////////////
    // Sampling Random Effects //
    /////////////////////////////
    uword n_samples = as<uword>(control["n_samples"]);
    uword n_iter = as<uword>(control["n_iter"]);
    uword n_b = b_mat.n_rows;
    mat scale_b = mat(n_b,  b_mat.n_cols, fill::ones) * 0.1;
    mat acceptance_b(n_b, b_mat.n_cols, fill::zeros);
    //
    field<vec> betas_it(betas.n_elem);
    mat out(n_b, n_samples);
    for (uword it = 0; it < n_samples; ++it) {
        vec bs_gammas_it = bs_gammas.col(it);
        vec gammas_it = gammas.col(it);
        vec alphas_it = alphas.col(it);
        for (uword i = 0; i < betas.n_elem; ++i) betas_it.at(i) = betas.at(i).col(it);
        vec sigmas_it = sigmas.col(it);
        mat L_it = L.slice(it);
        vec sds_it = sds.col(it);
        ///////////////////////
        vec W0H_bs_gammas = W0_H * bs_gammas_it;
        vec W0h_bs_gammas(W0_h.n_rows);
        vec W0H2_bs_gammas(W0_H2.n_rows);
        if (any_event) {
            W0h_bs_gammas = W0_h * bs_gammas_it;
        }
        if (any_interval) {
            W0H2_bs_gammas = W0_H2 * bs_gammas_it;
        }
        vec WH_gammas(W0_H.n_rows);
        vec Wh_gammas(W0_h.n_rows);
        vec WH2_gammas(W0_H2.n_rows);
        if (any_gammas) {
            WH_gammas = W_H * gammas_it;
        }
        if (any_gammas && any_event) {
            Wh_gammas = W_h * gammas_it;
        }
        if (any_gammas && any_interval) {
            WH2_gammas = W_H2 * gammas_it;
        }
        mat Wlong_H =
            calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, Wlong_sds, betas_it, b,
                            id_H_, FunForms, FunForms_ind, Funs_FunForms);
        vec WlongH_alphas = Wlong_H * alphas_it;
        mat Wlong_h(W0_h.n_rows, Wlong_H.n_cols);
        vec Wlongh_alphas(W0_h.n_rows);
        if (any_event) {
            Wlong_h =
                calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, Wlong_sds, betas_it, b,
                                id_h, FunForms, FunForms_ind, Funs_FunForms);
            Wlongh_alphas = Wlong_h * alphas_it;
        }
        mat Wlong_H2(W0_H2.n_rows, Wlong_H.n_cols);
        vec WlongH2_alphas(W0_H2.n_rows);
        if (any_interval) {
            Wlong_H2 =
                calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, Wlong_sds, betas_it,
                                b, id_H_, FunForms, FunForms_ind, Funs_FunForms);
            WlongH2_alphas = Wlong_H2 * alphas_it;
        }
        vec logLik_surv =
            log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                     WH_gammas, Wh_gammas, WH2_gammas,
                     WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                     log_Pwk, log_Pwk2, indFast_H, indFast_h,
                     which_event, which_right_event, which_left,
                     any_interval, which_interval);
        ///
        field<vec> eta_it = linpred_mixed(X, betas_it, Z, b, idL);
        vec logLik_long = log_long(y, eta_it, sigmas_it, extra_parms, families,
                                   links, ids, unq_idL, n_b);
        ///
        vec logLik_re = log_re(b_mat, L_it, sds_it);
        // calculate the denominator
        vec denominator_b =
            logLik_long + logLik_surv + logLik_re;
        for (uword j = 0; j < n_iter; ++j) {

        }
    }

    return(out);
}

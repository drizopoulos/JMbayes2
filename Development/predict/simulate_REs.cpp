#include <RcppArmadillo.h>
#include "JMbayes2_LogDens.h"
#include "JMbayes2_Funs.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
cube simulate_REs (List Data, List MCMC, List control) {
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
    //
    field<uvec> ind_RE = List2Field_uvec(as<List>(Data["ind_RE"]), true);
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
    mat bs_gammas = as<mat>(MCMC["bs_gammas"]);
    mat gammas = as<mat>(MCMC["gammas"]);
    mat alphas = as<mat>(MCMC["alphas"]);
    field<mat> betas = List2Field_mat(as<List>(MCMC["betas"]));
    mat sigmas = as<mat>(MCMC["sigmas"]);
    cube D = as<cube>(MCMC["D"]);

    //////////////////////////////
    // Sampling Random Effects //
    /////////////////////////////
    uword n_b = b_mat.n_rows;
    mat scale_b = mat(n_b,  b_mat.n_cols, fill::ones) * 0.1;
    mat acceptance_b(n_b, b_mat.n_cols, fill::zeros);

    cube out(1, 1, 1);
    return(D);
}

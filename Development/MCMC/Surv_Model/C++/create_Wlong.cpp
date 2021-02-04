#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

static double const Const_Unif_Proposal = 0.5 * std::pow(12.0, 0.5);
static double const log2pi = std::log(2.0 * M_PI);


field<mat> List2Field_mat (const List &Mats) {
    uword n_list = Mats.size();
    field<mat> res(n_list);
    for (uword i = 0; i < n_list; ++i) {
        res.at(i) = as<mat>(Mats[i]);
    }
    return res;
}

field<uvec> List2Field_uvec (const List &uVecs, bool substract1 = true) {
    uword n_list = uVecs.size();
    field<uvec> res(n_list);
    for (uword i = 0; i < n_list; ++i) {
        if (substract1) {
            res.at(i) = as<arma::uvec>(uVecs[i]) - 1;
        } else {
            res.at(i) = as<arma::uvec>(uVecs[i]);
        }
    }
    return res;
}

void transf_eta (mat &eta, const CharacterVector &fun_nams) {
    uword n = eta.n_cols;
    for (uword i = 0; i < n; i++) {
        if (fun_nams[i] == "identity") continue;
        if (fun_nams[i] == "expit") {
            eta.col(i) = 1.0 / (1.0 + trunc_exp(- eta.col(i)));
        } else if (fun_nams[i] == "exp" || fun_nams[i] == "dexp") {
            eta.col(i) = trunc_exp(eta.col(i));
        } else if (fun_nams[i] == "dexpit") {
            mat pp = 1.0 / (1.0 + trunc_exp(- eta.col(i)));
            eta.col(i) = pp % (1.0 - pp);
        }
    }
}

// [[Rcpp::export]]
field<mat> create_Wlong_cpp (List Data) {
    field<mat> eta = List2Field_mat(as<List>(Data["eta"]));
    field<mat> U = List2Field_mat(as<List>(Data["U"]));
    field<uvec> FunForms = List2Field_uvec(as<List>(Data["FunForms_cpp"]), true);
    field<uvec> ind = List2Field_uvec(as<List>(Data["FunForms_ind"]), true);
    List Funs_FunForms = as<List>(Data["Funs_FunForms"]);
    uword n_outcomes = eta.n_elem;
    field<mat> out(n_outcomes);
    for (uword i = 0; i < n_outcomes; ++i) {
        CharacterVector Funs_i = Funs_FunForms[i];
        mat eta_i = eta.at(i);
        transf_eta(eta_i, Funs_i);
        uvec FF_i = FunForms.at(i);
        mat U_i = U.at(i);
        uvec ind_i = ind.at(i);
        mat Wlong_i(eta_i.n_rows, U_i.n_cols, fill::ones);
        Wlong_i.cols(FF_i) %= eta_i.cols(ind_i);
        out.at(i) = U_i % Wlong_i;
    }
    return out;
}

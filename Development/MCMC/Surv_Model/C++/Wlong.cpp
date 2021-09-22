#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

field<mat> List2Field_mat (const List &Mats) {
    uword n_list = Mats.size();
    field<mat> res(n_list);
    for (uword i = 0; i < n_list; ++i) {
        res.at(i) = as<mat>(Mats[i]);
    }
    return res;
}

field<vec> List2Field_vec (const List &Vecs) {
    uword n_list = Vecs.size();
    field<vec> res(n_list);
    for (uword i = 0; i < n_list; ++i) {
        res.at(i) = as<vec>(Vecs[i]);
    }
    return res;
}

mat docall_cbindF (const field<mat> &Mats) {
    uword n = Mats.n_elem;
    uvec ncols(n);
    for (uword k = 0; k < n; ++k) {
        ncols.at(k) = Mats.at(k).n_cols;
    }
    uword N = sum(ncols);
    uword col_start = 0;
    uword col_end = ncols.at(0) - 1;
    mat out(Mats.at(0).n_rows, N);
    for (uword k = 0; k < n; ++k) {
        if (k > 0) {
            col_start += ncols.at(k - 1);
            col_end += ncols.at(k);
        }
        out.cols(col_start, col_end) = Mats.at(k);
    }
    return out;
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

// [[Rcpp::export]]
mat transf_etaC (const mat &eta, const CharacterVector &fun_nams) {
    uword k = fun_nams.length();
    mat out(eta.n_rows, k, fill::zeros);
    for (uword i = 0; i < k; i++) {
        if (fun_nams[i] == "identity") {
            out.col(i) = eta;
        } else if (fun_nams[i] == "expit") {
            out.col(i) = 1.0 / (1.0 + trunc_exp(- eta));
        } else if (fun_nams[i] == "exp" || fun_nams[i] == "dexp") {
            out.col(i) = trunc_exp(eta);
        } else if (fun_nams[i] == "dexpit") {
            mat pp = 1.0 / (1.0 + trunc_exp(- eta));
            out.col(i) = pp * (1.0 - pp);
        } else if (fun_nams[i] == "log") {
            out.col(i) = trunc_log(eta);
        } else if (fun_nams[i] == "log2") {
            out.col(i) = log2(eta);
        } else if (fun_nams[i] == "log10") {
            out.col(i) = log10(eta);
        } else if (fun_nams[i] == "sqrt") {
            out.col(i) = sqrt(eta);
        } else if (fun_nams[i] == "square") {
            out.col(i) = square(eta);
        } else if (fun_nams[i] == "cubic") {
            out.col(i) = square(eta) % eta;
        }
    }
    return out;
}

// [[Rcpp::export]]
field<mat> create_WlongC (const List &eta_, const List &U_,
                          const List &model_info) {
    field<mat> eta = List2Field_mat(eta_);
    field<mat> U = List2Field_mat(U_);
    field<uvec> FunForms = List2Field_uvec(as<List>(model_info["FunForms_cpp"]), true);
    field<uvec> FunForms_ind = List2Field_uvec(as<List>(model_info["FunForms_ind"]), true);
    List Funs_FunForms = as<List>(model_info["Funs_FunForms"]);
    //////////////////////////////////////
    uword n_outcomes = eta.n_elem;
    field<mat> out(n_outcomes);
    for (uword i = 0; i < n_outcomes; ++i) {
        mat eta_i = eta.at(i);
        List Funs_i = Funs_FunForms[i];
        uword n = Funs_i.length();
        field<mat> res(n);
        for (uword j = 0; j < n; ++j) {
            res.at(j) = transf_etaC(eta_i.col(j), Funs_i[j]);
        }
        mat Res = docall_cbindF(res);
        uvec FF_i = FunForms.at(i);
        uvec ind_i = FunForms_ind.at(i);
        mat U_i = U.at(i);
        mat Wlong_i(U_i.n_rows, U_i.n_cols, fill::ones);
        Wlong_i.cols(FF_i) %= Res;
        out.at(i) = U_i % Wlong_i;
    }
    return out;
}

// [[Rcpp::export]]
CharacterVector test (const List &model_info) {
    List Funs_FunForms = as<List>(model_info["Funs_FunForms"]);
    CharacterVector Funs_i = Funs_FunForms[1];
    return Funs_i;
}

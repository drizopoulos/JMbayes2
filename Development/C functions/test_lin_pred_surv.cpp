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
field<mat> linpred_surv2(List XX, List bbetas, List ZZ, List bb, const uvec &id) {
    field<mat> X = List2Field_mat(XX);
    field<mat> Z = List2Field_mat(ZZ);
    field<mat> b = List2Field_mat(bb);
    field<vec> betas = List2Field_vec(bbetas);
    uword n_outcomes = X.n_elem;
    field<mat> out(n_outcomes);
    for (uword i = 0; i < n_outcomes; ++i) {
        const mat& X_i = X.at(i);
        const vec& betas_i = betas.at(i);
        const mat& Z_i = Z.at(i);
        const mat& b_i = b.at(i);
        uword n_betas = betas_i.n_rows;
        uword n_REs = b_i.n_cols;
        uword n_forms = X_i.n_cols / n_betas;
        out.at(i).set_size(X_i.n_rows, n_forms);
        mat b_id = b_i.rows(id);
        for (uword j = 0; j < n_forms; ++j) {
            uword j_beta_start = j * n_betas;
            uword j_beta_end = j_beta_start + n_betas - 1;
            uword j_RE_start = j * n_REs;
            uword j_RE_end = j_RE_start + n_REs - 1;
            out.at(i).col(j) = X_i.cols(j_beta_start, j_beta_end) * betas_i +
                arma::sum(Z_i.cols(j_RE_start, j_RE_end) % b_id, 1);
        }
    }
    return out;
}

// [[Rcpp::export]]
field<mat> linpred_surv (List XX, List bbetas, List ZZ, List bb, const uvec &id) {
    field<mat> X = List2Field_mat(XX);
    field<mat> Z = List2Field_mat(ZZ);
    field<mat> b = List2Field_mat(bb);
    field<vec> betas = List2Field_vec(bbetas);
    uword n_outcomes = X.n_elem;
    field<mat> out(n_outcomes);
    for (uword i = 0; i < n_outcomes; ++i) {
        mat X_i = X.at(i);
        vec betas_i = betas.at(i);
        mat Z_i = Z.at(i);
        mat b_i = b.at(i);
        uword n_betas = betas_i.n_rows;
        uword n_REs = b_i.n_cols;
        uword n_forms = X_i.n_cols / n_betas;
        mat out_i(X_i.n_rows, n_forms);
        out.at(i) = out_i;
        for (uword j = 0; j < n_forms; ++j) {
            mat X_ij = X_i.cols(j * n_betas, (j + 1) * n_betas - 1);
            mat Z_ij = Z_i.cols(j * n_REs, (j + 1) * n_REs - 1);
            out.at(i).col(j) = X_ij * betas_i +
                arma::sum(Z_ij % b_i.rows(id), 1);
        }
    }
    return out;
}

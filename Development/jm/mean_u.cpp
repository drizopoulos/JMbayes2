#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
void mean_u_cpp (mat &mean_u_i, const vec &betas_i, const uvec &x_in_z_i,
                 const mat &Xbase_i, const uvec &baseline_i,
                 const uvec &unq_idL_i) {
    uword n = x_in_z_i.n_rows;
    if (mean_u_i.n_cols == n) {
        mean_u_i.each_row() = betas_i.elem(x_in_z_i).t();
    } else {
        mean_u_i.cols(0, n - 1).each_row() = betas_i.elem(x_in_z_i).t();
    }
    if (is_finite(Xbase_i)) {
        mean_u_i(unq_idL_i) = Xbase_i * betas_i.rows(baseline_i);
    }
}

// [[Rcpp::export]]
void update_mean_u (field<mat> mean_u, const field<vec> &betas,
                    const field<mat> &Xbase, const field<uvec> &x_in_z,
                    const field<uvec> &baseline, const field<uvec> &unq_idL) {
    uword n = mean_u.n_elem;
    for (uword i = 0; i < n; ++i) {
        vec betas_i = betas.at(i);
        mat Xbase_i = Xbase.at(i);
        uvec xinz_i = x_in_z.at(i);
        uvec base_i = baseline.at(i);
        uvec rowind_i = unq_idL.at(i);
        uword n = xinz_i.n_rows;
        if (mean_u.at(i).n_cols == n) {
            mean_u.at(i).each_row() = betas_i.elem(xinz_i).t();
        } else {
            mean_u.at(i).cols(0, n - 1).each_row() = betas_i.elem(xinz_i).t();
        }
        if (is_finite(Xbase_i)) {
            mean_u.at(i)(rowind_i) = Xbase_i * betas_i.rows(base_i);
        }
    }
}







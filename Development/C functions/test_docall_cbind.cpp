#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
mat docall_cbindF_fast(const field<mat> &Mats) {
    uword n = Mats.n_elem;
    uword total_cols = 0;
    for (uword k = 0; k < n; ++k) {
        total_cols += Mats.at(k).n_cols;
    }
    uword n_rows = Mats.at(0).n_rows;
    mat out(n_rows, total_cols);
    uword current_col = 0;
    for (uword k = 0; k < n; ++k) {
        const mat& M = Mats.at(k);
        uword m_cols = M.n_cols;
        if (m_cols > 0) {
            out.cols(current_col, current_col + m_cols - 1) = M;
            current_col += m_cols;
        }
    }
    return out;
}


// [[Rcpp::export]]
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



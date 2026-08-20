#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

static double const log2pi = std::log(2.0 * M_PI);

void inplace_UpperTrimat_mult (rowvec &x, const mat &trimat) {
    // in-place multiplication of x with an upper triangular matrix trimat
    // because in-place assignment is much faster but careful in use because
    // it changes the input vector x, i.e., not const
    uword const n = trimat.n_cols;
    for (uword j = n; j-- > 0;) {
        double tmp(0.0);
        for (uword i = 0; i <= j; ++i)
            tmp += trimat.at(i, j) * x.at(i);
        x.at(j) = tmp;
    }
}

// [[Rcpp::export]]
vec log_dmvnrm_chol (const mat &x, const mat &L) {
    // fast log density of the multivariate normal distribution
    // L is the Cholesky factor of the covariance matrix.
    using arma::uword;
    uword const n = x.n_rows, xdim = x.n_cols;
    vec out(n);
    mat V = inv(trimatu(L));
    double const log_det = sum(log(V.diag())),
        constants = -(double)xdim / 2.0 * log2pi,
        other_terms = constants + log_det;
    rowvec z_i(xdim);
    for (uword i = 0; i < n; i++) {
        z_i = x.row(i);
        inplace_UpperTrimat_mult(z_i, V);
        out.at(i) = other_terms - 0.5 * dot(z_i, z_i);
    }
    return out;
}

// [[Rcpp::export]]
vec log_dmvnrm_chol (const mat &x, const mat &L) {
    uword const n = x.n_rows, k = x.n_cols;
    vec out(n, arma::fill::none);
    mat V = inv(trimatu(L));
    double log_det = -arma::sum(arma::log(L.diag()));
    double constants = -(double)k / 2.0 * log2pi;
    double other_terms = constants + log_det;
    std::vector<double> z(k);
    for (uword i = 0; i < n; ++i) {
        for(uword j = 0; j < k; ++j) {
            z[j] = x.at(i, j);
        }
        for (uword j = k; j-- > 0;) {
            double tmp = 0.0;
            const double* col_j = V.colptr(j);
            for (uword c = 0; c <= j; ++c) {
                tmp += col_j[c] * z[c];
            }
            z[j] = tmp;
        }
        double sq_dist = 0.0;
        for(uword j = 0; j < k; ++j) {
            sq_dist += z[j] * z[j];
        }
        out[i] = other_terms - 0.5 * sq_dist;
    }
    return out;
}


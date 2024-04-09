#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

mat cov2cor (const mat &V) {
    vec Is = sqrt(1.0 / V.diag());
    mat out = V.each_col() % Is;
    out.each_row() %= Is.t();
    return out;
}

// [[Rcpp::export]]
mat LL(const mat &L, const uword &i, const umat &ind_zero_D) {
    uvec upper_part = trimatu_ind(size(L),  1);
    mat proposed_L(size(L), fill::zeros);
    vec l = L(upper_part);
    vec proposed_l = l;
    double ss = proposed_l(i) + 0.01;
    proposed_l(i) = ss;
    proposed_L(upper_part) = proposed_l;
    uword nn = ind_zero_D.n_rows;
    for (uword j = 0; j < nn; ++j) {
        uword j0 = ind_zero_D(j, 0);
        uword j1 = ind_zero_D(j, 1);
        mat copy_proposed_L = proposed_L;
        copy_proposed_L(j0, j1) = 0.0;
        proposed_L(j0, j1) = -sum(copy_proposed_L.col(j0) % copy_proposed_L.col(j1));
    }
    uword n = L.n_rows;
    for (uword j = 0; j < n; ++j) {
        vec ll = proposed_L.col(j);
        proposed_L(j, j) = sqrt(1 - dot(ll, ll));
    }
    return proposed_L;
}

arma::uvec std_setdiff(arma::uvec& x, arma::uvec& y) {
    std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
    std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
    std::vector<int> out;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(out, out.end()));
    return arma::conv_to<arma::uvec>::from(out);
}

// [[Rcpp::export]]
uvec ut_ind (const mat &L, const umat &ind_zero_D) {
    uvec upper_part = trimatu_ind(size(L),  1);
    uvec excl = ind_zero_D.col(0) + ind_zero_D.col(1) * L.n_rows;
    uvec out = std_setdiff(upper_part, excl);
    return out;
}


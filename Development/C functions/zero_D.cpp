#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

double logSumExp (const vec &x) {
    double maxval = max(x);
    double out = maxval + log(sum(exp(x - maxval)));
    return out;
}

// [[Rcpp::export]]
vec lseC (const vec &xx, const uvec &id_h2, const uvec &intgr_ind) {
    uvec unq_idh2 = unique(id_h2);
    uword n = unq_idh2.n_rows;
    vec out(n);
    for (uword i = 0; i < n; ++i) {
        uvec idx = find(id_h2 == i);
        uvec intgr_i = intgr_ind(idx);
        uword nn = intgr_i.n_rows;
        if (nn > 0) {
            out(i) = logSumExp(xx(idx));
        } else {
            out(i) = sum(xx(idx));
        }
    }
    vec res(xx.n_rows);
    res = out;
    return res;
}

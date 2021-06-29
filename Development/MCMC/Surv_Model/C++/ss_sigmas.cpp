#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


uword n_field (const field<vec> &x) {
    uword n = x.n_rows;
    uword out = 0;
    for (uword i = 0; i < n; ++i)
        out += x.at(i).n_rows;
    return out;
}

// [[Rcpp::export]]
field<vec> create_sigmas_field (const field<vec> &sigmas,
                                const uvec &ss_sigmas,
                                const field<uvec> &idL) {
    uword n = sigmas.n_rows;
    Rcpp::Rcout << "n is: " << n_field(sigmas);
    field<vec> out(n);
    for (uword i = 0; i < n; ++i) {
        vec sigmas_i = sigmas.at(i);
        uvec id_i = idL.at(i);
        if (ss_sigmas.at(i)) {
            out.at(i) = sigmas_i.rows(id_i);
        } else {
            vec xx(id_i.n_rows);
            xx.fill(as_scalar(sigmas_i));
            out.at(i) = xx;
        }
    }
    return out;
}




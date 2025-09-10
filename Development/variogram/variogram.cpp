#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
mat variogram_cpp(const field<vec> &y, const field<vec> &times) {
    uword n = y.n_elem;
    uvec Nis(n);
    for (uword i = 0; i < n; ++i) {
        uword n_i = y.at(i).n_rows;
        if (n_i > 1) {
            Nis.at(i) = Rf_choose(y.at(i).n_rows, 2);
        } else {
            Nis.at(i) = 1;
        }
    }
    mat out(sum(Nis), 2, fill::zeros);
    uword str = 0;
    uword end = Nis.at(0) - 1;
    for (uword i = 0; i < n; ++i) {
        vec y_i = y.at(i);
        vec times_i = times.at(i);
        uword n_i = y.at(i).n_rows;
        if (n_i > 1) {
            uword ncombs = Nis.at(i);
            vec diffs_i(ncombs, fill::zeros);
            vec lags_i(ncombs, fill::zeros);
            uword pos = 0;
            for (uword j = 0; j < n_i; ++j) {
                for (uword k = j + 1; k < n_i; ++k) {
                    double diff = y_i.at(j) - y_i.at(k);
                    diffs_i.at(pos) = 0.5 * diff * diff;
                    lags_i.at(pos) = abs(times_i.at(j) - times_i.at(k));
                    ++pos;
                }
            }
            out.submat(str, 0, end, 0) = lags_i;
            out.submat(str, 1, end, 1) = diffs_i;
        }
        if (i + 1 < n) {
            str += Nis.at(i);
            end += Nis.at(i + 1);
        }
    }
    return out;
}



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

// [[Rcpp::export]]
mat variogram_cpp_fast(const field<vec> &y, const field<vec> &times) {
    uword n = y.n_elem;

    // 1. FAST COMBINATORICS (Eliminate Rf_choose and the Nis vector)
    uword total_rows = 0;
    for (uword i = 0; i < n; ++i) {
        uword n_i = y.at(i).n_elem;
        if (n_i > 1) {
            total_rows += (n_i * (n_i - 1)) / 2; // Bare-metal integer math
        } else {
            total_rows += 1;
        }
    }

    mat out(total_rows, 2, fill::zeros);
    uword pos = 0; // Global position tracker

    for (uword i = 0; i < n; ++i) {

        // 2. REFERENCES TO PREVENT CLONING
        const vec& y_i = y.at(i);
        const vec& times_i = times.at(i);
        uword n_i = y_i.n_elem;

        if (n_i > 1) {
            for (uword j = 0; j < n_i; ++j) {

                // 3. INNER LOOP CACHING
                // Pull these out of the 'k' loop so the CPU registers them once
                double y_j = y_i.at(j);
                double t_j = times_i.at(j);

                for (uword k = j + 1; k < n_i; ++k) {
                    double diff = y_j - y_i.at(k);

                    // 4. DIRECT ASSIGNMENT (No temp vectors or submat copies)
                    out.at(pos, 0) = std::abs(t_j - times_i.at(k));
                    out.at(pos, 1) = 0.5 * diff * diff;

                    ++pos;
                }
            }
        } else {
            // Replicates your original logic: leaves a (0,0) row if n_i <= 1
            ++pos;
        }
    }

    return out;
}

// [[Rcpp::export]]
double total_var_cpp(const field<vec> &y) {
    uword n = y.n_elem;
    double out = 0.0;
    double s = 0.0;
    double count = 0.0;
    for (uword i = 0; i < n; ++i) {
        vec y_i = y.at(i);
        uword n_i = y.at(i).n_rows;
        for (uword j = i + 1; j < n; ++j) {
            vec y_j = y.at(j);
            uword n_j = y.at(j).n_rows;
            for (uword ii = 0; ii < n_i; ++ii) {
                for (uword jj = 0; jj < n_j; ++jj) {
                    double d = y_i.at(ii) - y_j.at(jj);
                    s += 0.5 * d * d;
                    count += 1;
                }
            }
        }
    }
    out = s / count;
    return out;
}



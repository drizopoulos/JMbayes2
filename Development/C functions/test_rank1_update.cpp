#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

static double const log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
mat rank1_update (const mat &U, // performs rank-1 update. If U = chol(M), returns chol(M + v * v.t())
                  const vec &v) {

    uword n = v.n_elem;
    mat Res = U;
    vec v2  = v;

    for (uword i = 0; i < n; i++) {
        double r = pow( pow(Res.at(i, i), 2) + pow(v2.at(i), 2), 0.5);
        double c = r / Res.at(i, i);
        double s = v2.at(i) / Res.at(i, i);
        Res.at(i, i) = r;

        if (i < n-1) {
            Res.submat(i, i + 1, i, n - 1) = (Res.submat(i, i + 1, i, n - 1) + s * v2.rows(i + 1, n - 1).t()) / c;
            v2.rows(i + 1, n - 1) = c * v2.rows(i + 1, n - 1) - s * Res.submat(i, i + 1, i, n - 1).t();
        }
    }
    return Res;
}

// [[Rcpp::export]]
mat rank1_update_fast(const mat &U, const vec &v) {
    uword n = v.n_elem;
    mat Res = U;
    vec v2 = v;
    for (uword i = 0; i < n; ++i) {
        double res_ii = Res.at(i, i);
        double v2_i   = v2.at(i);
        double r = std::sqrt(res_ii * res_ii + v2_i * v2_i);
        double c = r / res_ii;
        double s = v2_i / res_ii;
        Res.at(i, i) = r;
        for (uword j = i + 1; j < n; ++j) {
            double res_ij = Res.at(i, j);
            double v2_j   = v2.at(j);
            double new_res_ij = (res_ij + s * v2_j) / c;
            Res.at(i, j) = new_res_ij;
            v2.at(j) = c * v2_j - s * new_res_ij;
        }
    }
    return Res;
}

// [[Rcpp::export]]
mat chol_update(const mat &U, // If U = chol(M), returns chol(M.submat(keep, keep))
                const uvec &keep) { // keep must be a sorted vector, i.e, {2, 4, 5}, and counts from 0

    // to improve: later we can try to extend this approach further to obtain inv(U) from the required inv(U_i)

    uvec rem = regspace<uvec>(0,  U.n_cols - 1); rem.shed_rows(keep); // cols-rows to remove
    mat Res = U;
    uword n = rem.n_elem;

    for (uword i = 0; i < n; i++) { // rank-1 update for each col-row to be removed

        uword last_col = Res.n_cols - 1;

        if(rem.at(i) < last_col) {
            Res.submat(rem.at(i) + 1, rem.at(i) + 1, last_col, last_col) = rank1_update(Res.submat(rem.at(i) + 1, rem.at(i) + 1, last_col, last_col),
                       Res.submat(rem.at(i), rem.at(i) + 1, rem.at(i), last_col).t());
        }

        Res.shed_row(rem.at(i));
        Res.shed_col(rem.at(i));
        rem = rem - 1;
    }
    return Res;
}


// [[Rcpp::export]]
mat chol_update_fast (const mat &U, const uvec &keep) {
    uword N = U.n_cols;
    mat Res = U; // Only one deep copy made
    // 1. Generate the removal indices exactly once
    uvec rem = regspace<uvec>(0, N - 1);
    rem.shed_rows(keep);
    // 2. Process each removal index without ever shrinking the matrix
    for (uword i = 0; i < rem.n_elem; ++i) {
        uword k = rem.at(i);
        // 3. INLINED RANK-1 UPDATE
        // We map the submatrix directly to the global coordinates of Res.
        for (uword row_idx = k + 1; row_idx < N; ++row_idx) {
            double res_ii = Res.at(row_idx, row_idx);
            // Res.at(k, row_idx) perfectly acts as our 'v' vector!
            double v_i    = Res.at(k, row_idx);
            // Fast hypotenuse math (avoids pow)
            double r = std::sqrt(res_ii * res_ii + v_i * v_i);
            double c = r / res_ii;
            double s = v_i / res_ii;
            Res.at(row_idx, row_idx) = r;
            for (uword col_idx = row_idx + 1; col_idx < N; ++col_idx) {
                double res_ij = Res.at(row_idx, col_idx);
                double v_j    = Res.at(k, col_idx);
                double new_res_ij = (res_ij + s * v_j) / c;
                Res.at(row_idx, col_idx) = new_res_ij;
                // Overwrite the removed row to act as our updated 'v' for the next step
                Res.at(k, col_idx) = c * v_j - s * new_res_ij;
            }
        }
    }
    // 4. Extract the final submatrix in one single allocation
    return Res.submat(keep, keep);
}

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec group_sum2(const arma::vec& x, const arma::uvec& ind) {
    arma::uword m = ind.n_elem;
    // Allocate only the final output vector (size M)
    arma::vec out(m);
    if (m == 0) return out; // Safety check for empty index vector
    // Extract raw pointers for maximum access speed
    const double* p_x = x.memptr();
    const arma::uword* p_ind = ind.memptr();
    double* p_out = out.memptr();
    arma::uword start = 0;
    for (arma::uword i = 0; i < m; ++i) {
        arma::uword end = p_ind[i];
        // Safety bound check to prevent crashing R if an index is too large
        if (end >= x.n_elem) end = x.n_elem - 1;
        double current_sum = 0.0;
        // Sum the elements for the current segment
        for (arma::uword j = start; j <= end; ++j) {
            current_sum += p_x[j];
        }
        p_out[i] = current_sum;
        // Next group starts exactly one element after the current group ends
        start = end + 1;
    }
    return out;
}

// [[Rcpp::export]]
uvec create_fast_ind (const uvec &group) {
    uword l = group.n_rows;
    if (l == 1) return group - 1;
    uvec ind = find(group.rows(1, l - 1) != group.rows(0, l - 2));
    uword k = ind.n_rows;
    ind.insert_rows(k, 1);
    ind.at(k) = l - 1;
    return ind;
}

// [[Rcpp::export]]
vec group_sum(const vec& x, const uvec& ind) {
    vec cumsum_x = cumsum(x);
    vec out = cumsum_x.rows(ind);
    out.insert_rows(0, 1);
    out = diff(out);
    return out;
}

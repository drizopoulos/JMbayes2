#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

uvec create_fast_ind (const uvec &group) {
  uword l = group.n_rows;
  uvec ind = find(group.rows(1, l - 1) != group.rows(0, l - 2));
  uword k = ind.n_rows;
  ind.insert_rows(k, 1);
  ind.at(k) = l - 1;
  return ind;
}

// [[Rcpp::export]]
field<uvec> fooid (field<uvec> idL_lp) {
  field<uvec> idL_lp_fast(idL_lp.n_elem);
  for (uword i = 0; i < idL_lp.n_elem; ++i) {
    idL_lp_fast.at(i) = create_fast_ind(idL_lp.at(i));
  }
  return idL_lp_fast;
}


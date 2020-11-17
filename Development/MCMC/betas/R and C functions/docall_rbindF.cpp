#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec docall_rbindF (const field<vec> &F) { // binds a field of vectors into one vector
  uword n = F.n_elem;
  uword nrows = 0;
  uvec rows(n);
  for (uword i = 0; i < n; i++) {
    rows.at(i) = F.at(i).n_rows;
    nrows += rows.at(i);
  }
  vec V(nrows);
  uword ii = 0;
  for (uword i = 0; i < n; i++) {
    V.rows(ii, ii - 1 + rows.at(i)) = F.at(i);
    ii += rows.at(i);
  }
  return V;
}

/*** R

n_elem <- 10
list_vec <- lapply(seq_len(n_elem), function(n) rep(n, n))

all.equal(docall_rbindF(list_vec),
          matrix(unlist(list_vec), ncol=1))
*/

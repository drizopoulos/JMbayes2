#ifndef JMBAYES2CPP_H
#define JMBAYES2CPP_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

//////////////////////////////////////////////////////////////////////////////////////////////////
// Include Function/Class/Template declarations/definitions here to be used in multiple .cpp files
//////////////////////////////////////////////////////////////////////////////////////////////////

//rowsum_svft (OLD_VERSION)
arma::vec rowsum_by_group_fast (const vec& x, const uvec& group) {
  vec cum_x = cumsum(x);
  vec out = cum_x.elem(group);
  out.insert_rows(0, 1);
  unsigned int n = out.n_elem;
  out = out.rows(1, n - 1) - out.rows(0, n - 2);
  return(out);
}

// logsumexp versions of rowsum
// recursive version
arma::vec rowlogsumexp_recursive_by_group (const vec& x, const uvec& group) {
  // initiate out vector
  arma::vec out = x.elem(group);
  // number of unique groups 
  unsigned int k = group.n_elem;
  // create index version of group with one zero at start
  arma::uvec group_index = uvec(group.n_elem + 1).fill(0);
  // add correct elements to group index
  group_index.subvec(1, group_index.n_elem - 1) = group.subvec(0, group.n_elem - 1);
  // initiate field of vecs: one vec for each subject
  // and other variables for loop outside due to scoping
  arma::field<vec> X_i(k);
  arma::vec n_i(k);
  arma::vec Lk_i(k);
  arma::vec logsumexp_i(k);
  unsigned int i;
  // initiate loop over each group to apply recursive LSE 
  for (unsigned int i = 1; i <= k; ++i) {
    X_i(i - 1) = sort(x.rows(group_index(i - 1), group_index(i)), "descend");
    n_i(i - 1) = X_i(i - 1).n_elem;
    Lk_i(i - 1) = X_i(i - 1).at(0);
    for (unsigned int j = 0; j < n_i(i - 1) - 1; ++j) {
      logsumexp_i(i - 1) = std::max(X_i(i - 1)(j + 1), Lk_i(i - 1)) + log1p(exp(-abs(X_i(i - 1)(j + 1) - Lk_i(i - 1))));
    }
    out(i - 1) = logsumexp_i(i - 1);
  }
  return(out);  
}

#endif

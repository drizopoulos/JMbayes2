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
  vec out = x.elem(group);
  // number of unique groups 
  unsigned int k = group.n_elem;
  // create new vector for indexing in the loop
  vec group_index = group.insert_rows(0, 1);
  // initiate loop over each group to apply recursive LSE 
  for (int i = 1; i < k; ++i) {
    vec x_i = x.rows(group_index.elem(i - 1), group_index.elem(i));
    int n_i = x_i.n_elem;
    vec x_i_sorted = sort(x_i, 'descend');
    double Lk = x_i_sorted(0);
    for (int j = 0; j < n_i - 1; ++i) {
      double maxelem = xi_sorted(j + 1).insert_rows(1, Lk).max();
      double logelem = log1p(exp(-abs(xi_sorted(j + 1) - Lk)));
      double logsumexp_i = maxelem + logelem;
    }
    out(i - 1) = logsumexp_i;
  }
  return(out);
}

#endif

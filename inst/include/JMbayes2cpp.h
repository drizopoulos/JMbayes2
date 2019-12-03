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

////////////////////////////////
// logsumexp versions of rowsum
////////////////////////////////

// recursive version
arma::vec rowlogsumexp_recursive_by_group (const arma::vec& x, const arma::uvec& group) {
  // initiate out vector
  arma::vec out = x.elem(group);
  // number of unique groups 
  unsigned int k = group.n_elem;
  arma::uvec group_index = uvec(group.n_elem + 1).fill(0);
  // create new vector for indexing in the loop
  group_index.subvec(1, group_index.n_elem - 1) = group.subvec(0, group.n_elem - 1) + 1;
  arma::uvec group_index_start = group_index.head(group_index.n_elem - 1);
  arma::uvec group_index_end = group_index.tail(group_index.n_elem - 1) - 1;
  // initiate field of vecs for each subject
  arma::field<vec> X_i(k);
  arma::vec n_i(k);
  arma::vec Lk_i(k);
  unsigned int i;
  // initiate loop over each group to apply recursive LSE 
  // create field with loop
  for (unsigned int i = 1; i <= k; ++i) {
    X_i(i - 1) = sort(x.rows(group_index_start(i - 1), group_index_end(i - 1)), "descend");
    n_i(i - 1) = X_i(i - 1).n_elem;
    Lk_i(i - 1) = X_i(i - 1).at(0);
    for (unsigned int j = 0; j < n_i(i - 1) - 1; ++j) {
      Lk_i(i - 1) = std::max(X_i(i - 1)(j + 1), Lk_i(i - 1)) + log1p(exp(-std::abs(X_i(i - 1)(j + 1) - Lk_i(i - 1))));
    }
    out(i - 1) = Lk_i(i - 1);
  }
  return(out);  
}

////////////////////////////////////////////////////////////////////

////////////////
// R to cpp
////////////////

// linpred_mixed
arma::field<arma::vec> linpred_mixed(const arma::field<arma::mat>& X, const arma::field<arma::vec>& betas, 
                                     const arma::field<arma::mat>& Z, const arma::field<arma::mat>& b, 
                                     const arma::field<arma::uvec>& id) {
  int n_outcomes = X.n_elem;
  arma::field<arma::vec> out;
  for (int i = 0; i < n_outcomes; i++) {
    arma::mat X_i = X(i);
    arma::vec betas_i = betas(i);
    arma::mat Z_i = Z(i);
    arma::mat b_i = b(i);
    arma::uvec id_i = id(i);
    out(i) = X_i * betas_i + arma::sum(Z_i % b_i.rows(id_i), 1); 
  }
  return(out);
}

////////////////////////////////////////////////////////////////////////////////
// linpred_surv
// NOTES: 
// cURRENTLY WORKS WITH X_H & Z_H as list of arrays, where we have one array per outcome one slice per functional form
// NOTE THAT FOR IT TO WORK WE NEED TO ALSO CORRECT INDEXING VIA ID_H e.g lapply(id_H, FUN = function(x) x - 1)
arma::field<arma::mat> linpred_surv_array(const arma::field<arma::cube>& X, const arma::field<arma::vec>& betas, 
                                          const arma::field<arma::cube>& Z, const arma::field<arma::mat>& b, 
                                          const arma::field<arma::uvec>& id) {
  int n_outcomes = X.n_elem;
  arma::field<arma::mat> out(n_outcomes);
  for (int i = 0; i < n_outcomes; i++) {
    arma::cube X_i = X(i);
    arma::vec betas_i = betas(i);
    arma::cube Z_i = Z(i);
    arma::mat b_i = b(i);
    arma::uvec id_i = id(i);
    int n_forms = X_i.n_slices;
    arma::mat out_i(X_i.n_rows, n_forms);
    out(i) = out_i;
    for (int j = 0; j < n_forms; j++) {
      out(i).col(j) = X_i.slice(j) * betas_i + arma::sum(Z_i.slice(j) % b_i.rows(id_i), 1);
    }
  }
  return(out);
}


/////////////////////////////////////////////////////////////////////////////////////////////////
// Does not work. Needs to be within function scope
// OR find a way tyo export as global variables 
// family_map
// maps strings to integers for use of switch
std::map<std::string, int> family_map;
family_map.insert(std::make_pair('gaussian', 1));
family_map.insert(std::make_pair('binomial', 2));
family_map.insert(std::make_pair('poisson', 3));
family_map.insert(std::make_pair('negative.binomial', 4));
////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

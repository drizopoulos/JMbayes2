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


////////////////////////////////////////////////////////////////////////////////
// create_Wlong
// NOTES: 
// cURRENTLY WORKS WITH arma::uvec ns_functional_forms_per_outcome and arma::uvec maxs_functional_forms_per_outcome
arma::field<arma::mat> create_Wlong(arma::field<arma::mat> eta, 
                                    SEXP functional_forms_per_outcome, 
                                    arma::uvec ns_functional_forms_per_outcome, 
                                    arma::uvec maxs_functional_forms_per_outcome, 
                                    arma::field<arma::mat> U) { 
  int n_etas = eta.n_elem;
  List ffpo(functional_forms_per_outcome);
  List inner_ffpo(n_etas);
  arma::field<arma::mat> Wlong(n_etas);
  arma::mat eta_i;
  arma::mat U_i;
  arma::mat Wlong_i;
  arma::uvec ind;
  for (int i = 0; i < n_etas; i++) {
    eta_i = eta(i);
    U_i = U(i);
    Wlong_i = arma::ones<arma::mat>(eta_i.n_rows, maxs_functional_forms_per_outcome(i));
    for (int j = 0; j < ns_functional_forms_per_outcome(i); j++) {
      ind = as<arma::uvec>(as<List>(ffpo[i])[j]) - 1;
      Wlong_i.cols(ind) = Wlong_i.cols(ind) % repelem(eta_i.col(j), 1, ind.n_elem);
    }
    Wlong(i) = U_i % Wlong_i;
  }
  return(Wlong);
}

//////////////////////////////////////////////////////////
// calculate_u
//////////////////////////////////////////////////////////
arma::field<arma::mat> calculate_u(arma::field<arma::mat> Xhc, 
                                   arma::field<arma::uvec> columns_HC, 
                                   arma::field<arma::vec> betas, 
                                   arma::field<arma::mat> b, 
                                   arma::field<arma::uvec> unq_idL) {
  arma::field<arma::mat>u(b);
  int n = Xhc.n_elem;
  arma::mat Xhc_i;
  arma::uvec columns_HC_i;
  arma::vec betas_i;
  arma::mat b_i;
  arma::uvec unq_idL_i;
  arma::mat mean_b_i;
  int ncol_b_i;
  arma::uvec index;
  arma::uvec cindex;
  for (int i = 0; i < n; i++) {
    Xhc_i = Xhc(i);
    columns_HC_i = columns_HC(i);
    betas_i = betas(i);
    b_i = b(i);
    unq_idL_i = unq_idL(i);
    mean_b_i = b_i * 0;
    ncol_b_i = b_i.n_cols;
    for (int j = 0; j < ncol_b_i; j++) {
      index = find(columns_HC_i == j + 1);
      cindex = j;
      mean_b_i(unq_idL_i - 1, cindex) = Xhc_i.cols(index) * betas_i(index);
    }
    u(i) = b_i + mean_b_i;
  }
  return(u);
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



//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////// List 2 Field
//////////////////////////////////////////////////////////////////////////////////////////////////
arma::field<arma::vec> List2Field_vec (const List & Vecs) {
  int n_list = Vecs.size();
  arma::field<arma::vec> res(n_list);
  for (int i = 0; i < n_list; ++i) {
    res.at(i) = as<arma::vec>(Vecs[i]);
  }
  return(res);
}

arma::field<arma::uvec> List2Field_uvec (const List & uVecs, bool substract1 = true) {
  int n_list = uVecs.size();
  arma::field<arma::uvec> res(n_list);
  for (int i = 0; i < n_list; ++i) {
    if (substract1) {
      res.at(i) = as<arma::uvec>(uVecs[i]) - 1;
    } else {
      res.at(i) = as<arma::uvec>(uVecs[i]);
    }
  }
  return(res);
}

arma::field<arma::mat> List2Field_mat (const List & Mats) {
  int n_list = Mats.size();
  arma::field<arma::mat> res(n_list);
  for (int i = 0; i < n_list; ++i) {
    res.at(i) = as<arma::mat>(Mats[i]);
  }
  return(res);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif

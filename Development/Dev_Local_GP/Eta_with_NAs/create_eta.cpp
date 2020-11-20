//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

field<mat> List2Field_mat (const List &Mats) {
  uword n_list = Mats.size();
  field<mat> res(n_list);
  for (uword i = 0; i < n_list; ++i) {
    res.at(i) = as<mat>(Mats[i]);
  }
  return res;
}

field<vec> List2Field_vec (const List &Vecs) {
  uword n_list = Vecs.size();
  field<vec> res(n_list);
  for (uword i = 0; i < n_list; ++i) {
    res.at(i) = as<vec>(Vecs[i]);
  }
  return res;
}

field<uvec> List2Field_uvec (const List &uVecs, bool substract1 = true) {
  uword n_list = uVecs.size();
  field<uvec> res(n_list);
  for (uword i = 0; i < n_list; ++i) {
    if (substract1) {
      res.at(i) = as<arma::uvec>(uVecs[i]) - 1;
    } else {
      res.at(i) = as<arma::uvec>(uVecs[i]);
    }
  }
  return res;
}

// [[Rcpp::export]]
field<vec> linpred_mixed (List model_data, List model_info, List initial_values,
                          List priors, List control, List vcov_prop) {
  field<mat> X = List2Field_mat(as<List>(model_data["X"]));
  field<mat> Z = List2Field_mat(as<List>(model_data["Z"]));
  field<vec> betas = List2Field_vec(as<List>(initial_values["betas"]));
  field<mat> b = List2Field_mat(as<List>(initial_values["b"]));
  field<uvec> idL = List2Field_uvec(as<List>(model_data["idL"]), true);
  uword n_outcomes = X.n_elem;
  field<vec> out(n_outcomes);
  for (uword i = 0; i < n_outcomes; ++i) {
    mat X_i = X.at(i);
    vec betas_i = betas.at(i);
    mat Z_i = Z.at(i);
    mat b_i = b.at(i);
    uvec idL_i = idL.at(i);
    out.at(i) = X_i * betas_i + arma::sum(Z_i % b_i.rows(idL_i), 1);
  }
  return out;
}
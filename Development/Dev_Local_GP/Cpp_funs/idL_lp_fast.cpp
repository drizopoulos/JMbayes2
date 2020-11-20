//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

field<vec> List2Field_vec (const List &Vecs) {
  uword n_list = Vecs.size();
  field<vec> res(n_list);
  for (uword i = 0; i < n_list; ++i) {
    res.at(i) = as<vec>(Vecs[i]);
  }
  return res;
}

field<uvec> List2Field_uvec (const List & uVecs, bool substract1 = true) {
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

field<vec> foo (List etas, List etas_new, List id_ind) {
  Rcout << "OK1" << endl;
  field<vec> eta = List2Field_vec(etas);
  Rcout << "OK2" << endl;
  uword n = eta.n_elem;
  Rcout << "OK3" << endl;
  field<vec> eta_new = List2Field_vec(etas_new);
  Rcout << "OK4" << endl;
  field<uvec> id = List2Field_uvec(id_ind, 0);
  Rcout << "OK5" << endl;
  uword k = id.n_elem;
  for (uword j = 0; j < n; j++){
    for (uword i = 0; i < k; i++ ) {
      eta.at(j).elem(id.at(i)) = eta_new.at(j).elem(id.at(i));
    }
    //Rcout << eta.at(0).elem(id.at(0)) << endl;
  }
  return eta;
}


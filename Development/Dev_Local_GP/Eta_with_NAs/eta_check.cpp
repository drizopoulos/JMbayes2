//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


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

void eta_check(field<vec> &eta, const field<vec> &eta_prop, const uvec &which, 
              const List &List_of_List_id, const uword &n, List model_data, const List &unq_idL_lst) {
  field<uvec> idL = List2Field_uvec(as<List>(model_data["idL"]), true);
  field<uvec> idL_lst = List2Field_uvec(as<List>(unq_idL_lst), true);
  uword n_outcomes = eta.n_elem;
  //uvec eta_indx = unique(as<uvec>(as<List>(as<List>(List_of_List_id)[0])[1]));
  //Rcout << eta_indx.at(0) << endl;
  //Rcout << find(idL.at(0) == eta_indx.at(0)) << endl;
  //Rcout << find(idL.at(0) == 1) << endl;
  //vec id_eta;
  for (uword i = 0; i < n; i++) {
    if (which.at(i) == 1) {
      uvec out_indx = idL_lst.at(i);
      uword n_j = out_indx.n_elem;
      for (uword j = 0; j < n_j; j++) {
        //id_eta = as<vec>(as<List>(List_of_List_id)[j]);
        uword k = out_indx.at(j);
        uvec eta_indx = unique(as<uvec>(as<List>(as<List>(List_of_List_id)[k])[i]));
        //Rcout << "at" << eta_indx.at(0) << endl;
        //Rcout << "begin" << eta_indx.front() << endl;
        //Rcout << eta.at(0).rows(eta_indx) << endl;
        //eta.at(j).rows(find(idL.at(j) == as<uvec>(as<List>(as<List>(List_of_List_id)[j])[i]))) = eta_prop.at(j).rows(find(idL.at(j) == as<uvec>(as<List>(as<List>(List_of_List_id)[j])[i])));
        //Rcout << eta.at(j).rows(find(idL.at(j) == eta_indx)) << endl;
        eta.at(k).rows(find(idL.at(k) == eta_indx.front())) = eta_prop.at(k).rows(find(idL.at(k) == eta_indx.front()));
        //Rcout << as<uvec>(as<List>(as<List>(List_of_List_id)[0])[311]) << endl;
      }
    }
  }
}

// [[Rcpp::export]]
field<vec> eta_check_wrap(field<vec> &eta, const field<vec> &eta_prop, const uvec &which, 
                          const List &List_of_List_id, const uword &n, List model_data, const List &unq_idL_lst) {
  
  eta_check(eta, eta_prop, which, List_of_List_id, n, model_data, unq_idL_lst);
  return eta;
}
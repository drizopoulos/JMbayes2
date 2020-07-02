#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

arma::vec group_sum_cpp (const arma::vec& x, const arma::uvec& ind) {
    arma::vec cumsum_x = cumsum(x);
    arma::vec out = cumsum_x.elem(ind);
    out.insert_rows(0, 1);
    out = diff(out);
    return(out);
}

arma::field<arma::mat> List2Field_mat (const List& Mats) {
    int n_list = Mats.size();
    arma::field<arma::mat> res(n_list);
    for (int i = 0; i < n_list; ++i) {
        res.at(i) = as<arma::mat>(Mats[i]);
    }
    return(res);
}

arma::field<arma::vec> List2Field_vec (const List& Vecs) {
    int n_list = Vecs.size();
    arma::field<arma::vec> res(n_list);
    for (int i = 0; i < n_list; ++i) {
        res.at(i) = as<arma::vec>(Vecs[i]);
    }
    return(res);
}

arma::field<arma::mat> create_storage(const arma::field<arma::vec>& F,
                                      const int& n_iter) {
    int n = F.size();
    arma::field<arma::mat> out(n);
    for (int i = 0; i < n; ++i) {
        arma::vec aa = F.at(i);
        int n_i = aa.n_rows;
        arma::mat tt(n_iter, n_i, fill::zeros);
        out.at(i) = tt;
    }
    return(out);
}

arma::uvec create_fast_ind (const arma::uvec& group) {
    unsigned int l = group.n_rows;
    arma::uvec ind = find(group.rows(1, l - 1) != group.rows(0, l - 2));
    unsigned int k = ind.n_rows;
    ind.insert_rows(k, 1);
    ind.at(k) = l - 1;
    return(ind);
}

// [[Rcpp::export]]
double log_density_surv_cpp (const arma::vec& W0H_bs_gammas,
                             const arma::vec& WH_gammas,
                             const arma::vec& WlongH_alphas,
                             const arma::vec& W0h_bs_gammas,
                             const arma::vec& Wh_gammas,
                             const arma::vec& Wlongh_alphas,
                             const arma::vec& W0H2_bs_gammas,
                             const arma::vec& WH2_gammas,
                             const arma::vec& WlongH2_alphas,
                             const arma::vec& log_Pwk,
                             const arma::uvec& indFast_H,
                             const arma::uvec& which_event,
                             const arma::uvec& which_right_event) {
    arma::vec lambda_H = W0H_bs_gammas + WH_gammas + WlongH_alphas;
    arma::vec H = group_sum_cpp(exp(log_Pwk + lambda_H), indFast_H);
    int n = H.n_rows;
    arma::vec lambda_h(n);
    lambda_h.elem(which_event - 1) = W0h_bs_gammas.elem(which_event - 1) +
        Wh_gammas.elem(which_event - 1) + Wlongh_alphas.elem(which_event - 1);
    arma::vec log_Lik_surv(n);
    log_Lik_surv.elem(which_right_event - 1) = - H.elem(which_right_event - 1);
    log_Lik_surv.elem(which_event - 1) += lambda_h.elem(which_event - 1);
    return(sum(log_Lik_surv));
}

// [[Rcpp::export]]
List mcmc (List model_data, List model_info, List initial_values,
           List priors, List control) {
    // outcome vectors and design matrices
    arma::vec Time_right = as<vec>(model_data["Time_right"]);
    arma::vec Time_left = as<vec>(model_data["Time_left"]);
    arma::vec Time_start = as<vec>(model_data["Time_start"]);
    arma::vec delta = as<vec>(model_data["Time_start"]);
    arma::uvec which_event = as<uvec>(model_data["which_event"]) - 1;
    arma::uvec which_right = as<uvec>(model_data["which_right"]) - 1;
    arma::uvec which_right_event = join_cols(which_event, which_right);
    arma::uvec which_left = as<uvec>(model_data["which_left"]) - 1;
    arma::uvec which_interval = as<uvec>(model_data["which_interval"]) - 1;
    arma::mat W0_H = as<mat>(model_data["W0_H"]);
    arma::mat W0_h = as<mat>(model_data["W0_h"]);
    arma::mat W0_H2 = as<mat>(model_data["W0_H2"]);
    arma::mat W_H = as<mat>(model_data["W_H"]);
    arma::mat W_h = as<mat>(model_data["W_h"]);
    arma::mat W_H2 = as<mat>(model_data["W_H2"]);
    List Wlong_H_ = as<List>(model_data["Wlong_H"]);
    field<mat> Wlong_H = List2Field_mat(Wlong_H_);
    List Wlong_h_ = as<List>(model_data["Wlong_h"]);
    field<mat> Wlong_h = List2Field_mat(Wlong_h_);
    List Wlong_H2_ = as<List>(model_data["Wlong_H2"]);
    field<mat> Wlong_H2 = List2Field_mat(Wlong_H2_);
    // other information
    int n = as<int>(model_data["n"]);
    arma::uvec idT = as<uvec>(model_data["idT"]) - 1;
    arma::vec log_Pwk = as<vec>(model_data["log_Pwk"]);
    arma::vec log_Pwk2 = as<vec>(model_data["log_Pwk2"]);
    arma::uvec id_H = create_fast_ind(as<uvec>(model_data["id_H"]));
    arma::uvec id_H2 = create_fast_ind(as<uvec>(model_data["id_H2"]));
    // initial values
    List betas_ = as<List>(initial_values["betas"]);
    field<vec> betas = List2Field_vec(betas_);
    List b_ = as<List>(initial_values["b"]);
    field<mat> b = List2Field_mat(b_);
    arma::vec bs_gammas = as<vec>(initial_values["bs_gammas"]);
    arma::vec gammas = as<vec>(initial_values["gammas"]);
    List alphas_ = as<List>(initial_values["alphas"]);
    field<vec> alphas = List2Field_vec(alphas_);
    // MCMC settings
    int n_iter = as<int>(control["n_iter"]);
    int n_burnin = as<int>(control["n_burnin"]);
    // scales
    arma::vec scale_bs_gammas(bs_gammas.n_rows, fill::ones);
    scale_bs_gammas = 0.1 * scale_bs_gammas;
    arma::vec scale_gammas(gammas.n_rows, fill::ones);
    scale_gammas = 0.1 * scale_gammas;
    // store results
    arma::mat res_bs_gammas(n_iter, bs_gammas.n_rows);
    arma::mat acceptance_bs_gammas(n_iter, bs_gammas.n_rows);
    arma::mat res_gammas(n_iter, gammas.n_rows);
    arma::mat acceptance_gammas(n_iter, gammas.n_rows);
    arma::field<mat> res_alphas = create_storage(alphas, n_iter);
    arma::field<mat> acceptance_alphas = create_storage(alphas, n_iter);
    res_alphas.at(0).at(0, 0) = 10;
    res_alphas.at(0).at(1, 0) = 20;
    res_alphas.at(1).at(0, 0) = 100;
    res_alphas.at(1).at(0, 1) = 200;
    res_alphas.at(1).at(0, 2) = 300;
    res_alphas.at(1).at(1, 0) = 101;
    res_alphas.at(1).at(1, 1) = 201;
    res_alphas.at(1).at(1, 2) = 301;


    return List::create(
        Named("Time_right") = Time_right,
        Named("Time_left") = Time_left,
        Named("delta") = delta,
        Named("which_event") = which_event + 1,
        Named("id_H") = res_alphas
    );

}




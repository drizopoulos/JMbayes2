#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec group_sum_cpp (const arma::vec& x, const arma::uvec& ind) {
    arma::vec cumsum_x = cumsum(x);
    arma::vec out = cumsum_x.elem(ind);
    out.insert_rows(0, 1);
    out = diff(out);
    return(out);
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

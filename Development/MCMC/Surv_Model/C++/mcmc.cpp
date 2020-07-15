#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

arma::vec group_sum (const arma::vec& x, const arma::uvec& ind) {
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

arma::vec Wlong_alphas_fun (const arma::field<arma::mat>& Mats,
                            const arma::field<arma::vec>& coefs) {
    int n = Mats.size();
    int n_rows = Mats.at(0).n_rows;
    arma::vec out(n_rows, fill::zeros);
    for (int k = 0; k < n; ++k) {
        out += Mats.at(k) * coefs.at(k);
    }
    return(out);
}

arma::mat docall_cbind (const List& Mats_) {
    field<mat> Mats = List2Field_mat(Mats_);
    int n = Mats.size();
    arma::uvec ncols(n);
    for (int k = 0; k < n; ++k) {
        ncols.at(k) = Mats.at(k).n_cols;
    }
    int N = sum(ncols);
    int col_start = 0;
    int col_end = ncols.at(0) - 1;
    arma::mat out(Mats.at(0).n_rows, N);
    for (int k = 0; k < n; ++k) {
        if (k > 0) {
            col_start += ncols.at(k - 1);
            col_end += ncols.at(k);
        }
        out.cols(col_start, col_end) = Mats.at(k);
    }
    return out;
}

arma::uvec create_fast_ind (const arma::uvec& group) {
    unsigned int l = group.n_rows;
    arma::uvec ind = find(group.rows(1, l - 1) != group.rows(0, l - 2));
    unsigned int k = ind.n_rows;
    ind.insert_rows(k, 1);
    ind.at(k) = l - 1;
    return(ind);
}

double logPrior(const vec& x, const vec& mean, const mat& Tau, double tau = 1) {
    vec z = (x - mean);
    return(- 0.5 * tau * as_scalar(z.t() * Tau * z));
}

arma::vec propose (const arma::vec& thetas, const arma::vec& scale, const int& i) {
    arma::vec proposed_thetas = thetas;
    proposed_thetas.at(i) = R::rnorm(thetas.at(i), scale.at(i));
    return proposed_thetas;
}

arma::field<arma::vec> propose_field (const arma::field<arma::vec>& thetas,
                                      const arma::field<arma::vec>& scale,
                                      const int& k, const int& i) {
    arma::field<arma::vec> proposed_thetas = thetas;
    proposed_thetas.at(k).at(i) = R::rnorm(thetas.at(k).at(i),
                                           scale.at(k).at(i));
    return proposed_thetas;
}

double robbins_monro (const double& scale, const double& acceptance_it,
                      const int& it, const double& target_acceptance = 0.45) {
    double step_length = scale / (target_acceptance * (1 - target_acceptance));
    double out;
    if (acceptance_it > 0) {
        out = scale + step_length * (1 - target_acceptance) / (double)it;
    } else {
        out = scale - step_length * target_acceptance / (double)it;
    }
    return out;
}

// [[Rcpp::export]]
double log_density_surv (const arma::vec& W0H_bs_gammas,
                         const arma::vec& W0h_bs_gammas,
                         const arma::vec& W0H2_bs_gammas,
                         const arma::vec& WH_gammas,
                         const arma::vec& Wh_gammas,
                         const arma::vec& WH2_gammas,
                         const arma::vec& WlongH_alphas,
                         const arma::vec& Wlongh_alphas,
                         const arma::vec& WlongH2_alphas,
                         const arma::vec& log_Pwk, const arma::vec& log_Pwk2,
                         const arma::uvec& indFast_H,
                         const arma::uvec& which_event,
                         const arma::uvec& which_right_event,
                         const arma::uvec& which_left,
                         const bool& any_interval,
                         const arma::uvec& which_interval) {
    arma::vec lambda_H = W0H_bs_gammas + WH_gammas + WlongH_alphas;
    arma::vec H = group_sum(exp(log_Pwk + lambda_H), indFast_H);
    int n = H.n_rows;
    arma::vec lambda_h(n);
    lambda_h.elem(which_event) = W0h_bs_gammas.elem(which_event) +
        Wh_gammas.elem(which_event) + Wlongh_alphas.elem(which_event);
    arma::vec log_Lik_surv(n);
    log_Lik_surv.elem(which_right_event) = - H.elem(which_right_event);
    log_Lik_surv.elem(which_event) += lambda_h.elem(which_event);
    log_Lik_surv.elem(which_left) = log1p(- exp(- H.elem(which_left)));
    arma::vec lambda_H2(lambda_H.n_rows);
    arma::vec H2(n);
    if (any_interval) {
        lambda_H2 = W0H2_bs_gammas + WH2_gammas + WlongH2_alphas;
        H2 = group_sum(exp(log_Pwk2 + lambda_H2), indFast_H);
        log_Lik_surv.elem(which_interval) = - H.elem(which_interval) +
            log(- expm1(- H2.elem(which_interval)));
    }
    return(sum(log_Lik_surv));
}

void update_bs_gammas (vec& bs_gammas, vec& gammas, vec& alphas,
                       vec& W0H_bs_gammas, vec& W0h_bs_gammas, vec& W0H2_bs_gammas,
                       vec& WH_gammas, vec& Wh_gammas, vec& WH2_gammas,
                       vec& WlongH_alphas, vec& Wlongh_alphas, vec& WlongH2_alphas,
                       vec& log_Pwk, vec& log_Pwk2, uvec& id_H,
                       uvec& which_event, uvec& which_right_event, uvec& which_left, uvec& which_interval,
                       bool& any_event, bool& any_interval,
                       vec& prior_mean_bs_gammas, mat& prior_Tau_bs_gammas, double& tau_bs_gammas,
                       vec& prior_mean_gammas, mat& prior_Tau_gammas,
                       double& denominator_surv, int& it,
                       /////
                       mat& W0_H, mat& W0_h, mat& W0_H2,
                       vec& scale_bs_gammas,
                       mat& acceptance_bs_gammas,
                       mat& res_bs_gammas) {
    for (unsigned int i = 0; i < bs_gammas.n_rows; ++i) {
        vec proposed_bs_gammas = propose(bs_gammas, scale_bs_gammas, i);
        vec proposed_W0H_bs_gammas = W0_H * proposed_bs_gammas;
        vec proposed_W0h_bs_gammas(W0_h.n_rows);
        vec proposed_W0H2_bs_gammas(W0_H2.n_rows);
        if (any_event) {
            proposed_W0h_bs_gammas = W0_h * proposed_bs_gammas;
        }
        if (any_interval) {
            proposed_W0H2_bs_gammas = W0_H2 * proposed_bs_gammas;
        }
        double numerator_surv =
            log_density_surv(proposed_W0H_bs_gammas, proposed_W0h_bs_gammas, proposed_W0H2_bs_gammas,
                             WH_gammas, Wh_gammas, WH2_gammas,
                             WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                             log_Pwk, log_Pwk2, id_H,
                             which_event, which_right_event, which_left,
                             any_interval, which_interval) +
                                 logPrior(proposed_bs_gammas, prior_mean_bs_gammas,
                                          prior_Tau_bs_gammas, tau_bs_gammas) +
                                              logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1);
        double log_ratio = numerator_surv - denominator_surv;
        if (is_finite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
            bs_gammas = proposed_bs_gammas;
            W0H_bs_gammas = proposed_W0H_bs_gammas;
            if (any_event) {
                W0h_bs_gammas = proposed_W0h_bs_gammas;
            }
            if (any_interval) {
                W0H2_bs_gammas = proposed_W0H2_bs_gammas;
            }
            denominator_surv = numerator_surv;
            acceptance_bs_gammas.at(it, i) = 1;
        }
        if (it > 19) {
            scale_bs_gammas.at(i) =
                robbins_monro(scale_bs_gammas.at(i),
                              acceptance_bs_gammas.at(it, i), it);
        }
        res_bs_gammas.at(it, i) = bs_gammas.at(i);
    }
}

void update_gammas (vec& bs_gammas, vec& gammas, vec& alphas,
                    vec& W0H_bs_gammas, vec& W0h_bs_gammas, vec& W0H2_bs_gammas,
                    vec& WH_gammas, vec& Wh_gammas, vec& WH2_gammas,
                    vec& WlongH_alphas, vec& Wlongh_alphas, vec& WlongH2_alphas,
                    vec& log_Pwk, vec& log_Pwk2, uvec& id_H,
                    uvec& which_event, uvec& which_right_event, uvec& which_left, uvec& which_interval,
                    bool& any_event, bool& any_interval,
                    vec& prior_mean_bs_gammas, mat& prior_Tau_bs_gammas, double& tau_bs_gammas,
                    vec& prior_mean_gammas, mat& prior_Tau_gammas,
                    double& denominator_surv, int& it,
                    /////
                    mat& W_H, mat& W_h, mat& W_H2,
                    vec& scale_gammas,
                    mat& acceptance_gammas,
                    mat& res_gammas) {
    for (unsigned int i = 0; i < gammas.n_rows; ++i) {
        vec proposed_gammas = propose(gammas, scale_gammas, i);
        vec proposed_WH_gammas = W_H * proposed_gammas;
        vec proposed_Wh_gammas(W_h.n_rows);
        vec proposed_WH2_gammas(W_H2.n_rows);
        if (any_event) {
            proposed_Wh_gammas = W_h * proposed_gammas;
        }
        if (any_interval) {
            proposed_WH2_gammas = W_H2 * proposed_gammas;
        }
        double numerator_surv =
            log_density_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                             proposed_WH_gammas, proposed_Wh_gammas, proposed_WH2_gammas,
                             WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                             log_Pwk, log_Pwk2, id_H,
                             which_event, which_right_event, which_left,
                             any_interval, which_interval) +
            logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas,
                     tau_bs_gammas) +
            logPrior(proposed_gammas, prior_mean_gammas, prior_Tau_gammas, 1);
        double log_ratio = numerator_surv - denominator_surv;
        if (is_finite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
            gammas = proposed_gammas;
            WH_gammas = proposed_WH_gammas;
            if (any_event) {
                Wh_gammas = proposed_Wh_gammas;
            }
            if (any_interval) {
                WH2_gammas = proposed_WH2_gammas;
            }
            denominator_surv = numerator_surv;
            acceptance_gammas.at(it, i) = 1;
        }
        if (it > 19) {
            scale_gammas.at(i) =
                robbins_monro(scale_gammas.at(i),
                              acceptance_gammas.at(it, i), it);
        }
        // store results
        res_gammas.at(it, i) = gammas.at(i);
    }
}


void update_alphas (vec& bs_gammas, vec& gammas, vec& alphas,
                    vec& W0H_bs_gammas, vec& W0h_bs_gammas, vec& W0H2_bs_gammas,
                    vec& WH_gammas, vec& Wh_gammas, vec& WH2_gammas,
                    vec& WlongH_alphas, vec& Wlongh_alphas, vec& WlongH2_alphas,
                    vec& log_Pwk, vec& log_Pwk2, uvec& id_H,
                    uvec& which_event, uvec& which_right_event, uvec& which_left, uvec& which_interval,
                    bool& any_event, bool& any_interval,
                    vec& prior_mean_bs_gammas, mat& prior_Tau_bs_gammas, double& tau_bs_gammas,
                    vec& prior_mean_gammas, mat& prior_Tau_gammas,
                    double& denominator_surv, int& it,
                    /////
                    mat& Wlong_H, mat& Wlong_h, mat& Wlong_H2,
                    vec& scale_alphas,
                    mat& acceptance_alphas,
                    mat& res_alphas) {
    for (unsigned int i = 0; i < alphas.n_rows; ++i) {
        vec proposed_alphas = propose(alphas, scale_alphas, i);
        vec proposed_WlongH_alphas = Wlong_H * proposed_alphas;
        vec proposed_Wlongh_alphas(Wlong_h.n_rows);
        if (any_event) {
            proposed_Wlongh_alphas = Wlong_h * proposed_alphas;
        }
        vec proposed_WlongH2_alphas(Wlong_H2.n_rows);
        if (any_interval) {
            proposed_WlongH2_alphas = Wlong_H2 * proposed_alphas;
        }
        double numerator_surv =
            log_density_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                             WH_gammas, Wh_gammas, WH2_gammas,
                             proposed_WlongH_alphas, proposed_Wlongh_alphas, proposed_WlongH2_alphas,
                             log_Pwk, log_Pwk2, id_H,
                             which_event, which_right_event, which_left,
                             any_interval, which_interval) +
                                 logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
                                 logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1);
        double log_ratio = numerator_surv - denominator_surv;
        if (is_finite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
            alphas = proposed_alphas;
            WlongH_alphas = proposed_WlongH_alphas;
            if (any_event) {
                Wlongh_alphas = proposed_Wlongh_alphas;
            }
            if (any_interval) {
                WlongH2_alphas = proposed_WlongH2_alphas;
            }
            denominator_surv = numerator_surv;
            acceptance_alphas.at(it, i) = 1;
        }
        if (it > 19) {
            scale_alphas.at(i) =
                robbins_monro(scale_alphas.at(i),
                              acceptance_alphas.at(it, i), it);
        }
        // store results
        res_alphas.at(it, i) = alphas.at(i);
    }
}

// [[Rcpp::export]]
List mcmc (List model_data, List model_info, List initial_values,
           List priors, List control) {
    // outcome vectors and design matrices
    vec Time_right = as<vec>(model_data["Time_right"]);
    vec Time_left = as<vec>(model_data["Time_left"]);
    vec Time_start = as<vec>(model_data["Time_start"]);
    vec delta = as<vec>(model_data["Time_start"]);
    uvec which_event = as<uvec>(model_data["which_event"]) - 1;
    uvec which_right = as<uvec>(model_data["which_right"]) - 1;
    uvec which_right_event = join_cols(which_event, which_right);
    uvec which_left = as<uvec>(model_data["which_left"]) - 1;
    uvec which_interval = as<uvec>(model_data["which_interval"]) - 1;
    mat W0_H = as<mat>(model_data["W0_H"]);
    mat W0_h = as<mat>(model_data["W0_h"]);
    mat W0_H2 = as<mat>(model_data["W0_H2"]);
    mat W_H = as<mat>(model_data["W_H"]);
    mat W_h = as<mat>(model_data["W_h"]);
    mat W_H2 = as<mat>(model_data["W_H2"]);
    mat W_bar = as<mat>(model_data["W_bar"]);
    List Wlong_H_ = as<List>(model_data["Wlong_H"]);
    mat Wlong_H = docall_cbind(Wlong_H_);
    List Wlong_h_ = as<List>(model_data["Wlong_h"]);
    mat Wlong_h = docall_cbind(Wlong_h_);
    List Wlong_H2_ = as<List>(model_data["Wlong_H2"]);
    mat Wlong_H2 = docall_cbind(Wlong_H2_);
    List Wlong_bar_ = as<List>(model_data["Wlong_bar"]);
    field<mat> Wlong_bar = List2Field_mat(Wlong_bar_);
    // other information
    uvec idT = as<uvec>(model_data["idT"]) - 1;
    vec log_Pwk = as<vec>(model_data["log_Pwk"]);
    vec log_Pwk2 = as<vec>(model_data["log_Pwk2"]);
    uvec id_H = create_fast_ind(as<uvec>(model_data["id_H"]));
    uvec id_H2 = create_fast_ind(as<uvec>(model_data["id_H2"]));
    bool any_gammas = as<bool>(model_data["any_gammas"]);
    bool any_event = which_event.n_rows > 0;
    bool any_interval = which_interval.n_rows > 0;
    // initial values
    List betas_ = as<List>(initial_values["betas"]);
    field<vec> betas = List2Field_vec(betas_);
    List b_ = as<List>(initial_values["b"]);
    field<mat> b = List2Field_mat(b_);
    vec bs_gammas = as<vec>(initial_values["bs_gammas"]);
    vec gammas = as<vec>(initial_values["gammas"]);
    vec alphas = as<vec>(initial_values["alphas"]);
    double tau_bs_gammas = as<double>(initial_values["tau_bs_gammas"]);
    // MCMC settings
    int n_iter = as<int>(control["n_iter"]);
    int n_burnin = as<int>(control["n_burnin"]);
    // priors
    vec prior_mean_bs_gammas = as<vec>(priors["mean_bs_gammas"]);
    mat prior_Tau_bs_gammas = as<mat>(priors["Tau_bs_gammas"]);
    vec prior_mean_gammas = as<vec>(priors["mean_gammas"]);
    mat prior_Tau_gammas = as<mat>(priors["Tau_gammas"]);
    vec prior_mean_alphas = as<vec>(priors["mean_alphas"]);
    mat prior_Tau_alphas = as<mat>(priors["Tau_alphas"]);
    double post_A_tau_bs_gammas = as<double>(priors["A_tau_bs_gammas"]) +
        0.5 * as<double>(priors["rank_Tau_bs_gammas"]);
    double prior_B_tau_bs_gammas = as<double>(priors["B_tau_bs_gammas"]);
    // scales
    vec scale_bs_gammas(bs_gammas.n_rows, fill::ones);
    scale_bs_gammas = 0.1 * scale_bs_gammas;
    vec scale_gammas(gammas.n_rows, fill::ones);
    scale_gammas = 0.1 * scale_gammas;
    vec scale_alphas(alphas.n_rows, fill::ones);
    scale_alphas = 0.1 * scale_alphas;
    // store results
    int n_bs_gammas = bs_gammas.n_rows;
    int n_gammas = gammas.n_rows;
    int n_alphas = alphas.n_rows;
    mat res_bs_gammas(n_iter, n_bs_gammas);
    mat acceptance_bs_gammas(n_iter, n_bs_gammas, fill::zeros);
    mat res_gammas(n_iter, n_gammas);
    vec res_W_bar_gammas(n_iter);
    mat acceptance_gammas(n_iter, n_gammas, fill::zeros);
    mat res_alphas(n_iter, n_alphas);
    mat acceptance_alphas(n_iter, n_alphas, fill::zeros);
    vec res_tau_bs_gammas(n_iter);
    // preliminaries
    vec W0H_bs_gammas = W0_H * bs_gammas;
    vec W0h_bs_gammas(W0_h.n_rows);
    vec W0H2_bs_gammas(W0_H2.n_rows);
    if (any_event) {
        W0h_bs_gammas = W0_h * bs_gammas;
    }
    if (any_interval) {
        W0H2_bs_gammas = W0_H2 * bs_gammas;
    }
    vec WH_gammas(W0_H.n_rows);
    vec Wh_gammas(W0_h.n_rows);
    vec WH2_gammas(W0_H2.n_rows);
    if (any_gammas) {
        WH_gammas = W_H * gammas;
    }
    if (any_gammas && any_event) {
        Wh_gammas = W_h * gammas;
    }
    if (any_gammas && any_interval) {
        WH2_gammas = WH2_gammas * gammas;
    }
    vec WlongH_alphas = Wlong_H * alphas;
    vec Wlongh_alphas(W0_h.n_rows);
    if (any_event) {
        Wlongh_alphas = Wlong_h * alphas;
    }
    vec WlongH2_alphas(W0_H2.n_rows);
    if (any_interval) {
        WlongH2_alphas = Wlong_H2 * alphas;
    }
    double denominator_surv =
        log_density_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                         WH_gammas, Wh_gammas, WH2_gammas,
                         WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                         log_Pwk, log_Pwk2, id_H,
                         which_event, which_right_event, which_left,
                         any_interval, which_interval) +
        logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
        logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1);
    for (int it = 0; it < n_iter; ++it) {
        update_bs_gammas(bs_gammas, gammas, alphas,
            W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
            WH_gammas, Wh_gammas, WH2_gammas,
            WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
            log_Pwk, log_Pwk2, id_H,
            which_event, which_right_event, which_left, which_interval,
            any_event, any_interval,
            prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas,
            prior_mean_gammas, prior_Tau_gammas,
            denominator_surv, it,
            /////
            W0_H, W0_h, W0_H2, scale_bs_gammas, acceptance_bs_gammas,
            res_bs_gammas);
        ////////////////////////////////////////////////////////////////////////
        double post_B_tau_bs_gammas = prior_B_tau_bs_gammas +
            0.5 * as_scalar(bs_gammas.t() * prior_Tau_bs_gammas * bs_gammas);
        tau_bs_gammas = R::rgamma(post_A_tau_bs_gammas, 1 / post_B_tau_bs_gammas);
        res_tau_bs_gammas.at(it) = tau_bs_gammas;
        ////////////////////////////////////////////////////////////////////////
        if (any_gammas) {
            update_gammas(bs_gammas, gammas, alphas,
                          W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                          WH_gammas, Wh_gammas, WH2_gammas,
                          WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                          log_Pwk, log_Pwk2, id_H,
                          which_event, which_right_event, which_left, which_interval,
                          any_event, any_interval,
                          prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas,
                          prior_mean_gammas, prior_Tau_gammas,
                          denominator_surv, it,
                          /////
                          W_H, W_h, W_H2, scale_gammas, acceptance_gammas,
                          res_gammas);
            res_W_bar_gammas.at(it) = as_scalar(W_bar * gammas);
        }
        ////////////////////////////////////////////////////////////////////////
        update_alphas(bs_gammas, gammas, alphas,
                      W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                      WH_gammas, Wh_gammas, WH2_gammas,
                      WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                      log_Pwk, log_Pwk2, id_H,
                      which_event, which_right_event, which_left, which_interval,
                      any_event, any_interval,
                      prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas,
                      prior_mean_gammas, prior_Tau_gammas,
                      denominator_surv, it,
                      /////
                      Wlong_H, Wlong_h, Wlong_H2, scale_alphas,
                      acceptance_alphas, res_alphas);
    }
    return List::create(
        Named("mcmc") = List::create(
            Named("bs_gammas") = res_bs_gammas.rows(n_burnin, n_iter - 1),
            Named("tau_bs_gammas") = res_tau_bs_gammas.rows(n_burnin, n_iter - 1),
            Named("gammas") = res_gammas.rows(n_burnin, n_iter - 1),
            Named("W_bar_gammas") = res_W_bar_gammas.rows(n_burnin, n_iter - 1),
            Named("alphas") = res_alphas.rows(n_burnin, n_iter - 1)
        ),
        Named("acc_rate") = List::create(
            Named("bs_gammas") = acceptance_bs_gammas.rows(n_burnin, n_iter - 1),
            Named("gammas") = acceptance_gammas.rows(n_burnin, n_iter - 1),
            Named("alphas") = acceptance_alphas.rows(n_burnin, n_iter - 1)
        )
    );
}




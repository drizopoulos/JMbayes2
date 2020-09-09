#include <RcppArmadillo.h>
#include "JMbayes2_D.h"
#include "JMbayes2_Surv.h"
# include "JMbayes2_Long.h"

// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List mcmc_cpp (List model_data, List model_info, List initial_values,
               List priors, List control, List vcov_prop) {
  // outcome vectors and design matrices
  vec Time_right = as<vec>(model_data["Time_right"]);
  vec Time_left = as<vec>(model_data["Time_left"]);
  vec Time_start = as<vec>(model_data["Time_start"]);
  vec delta = as<vec>(model_data["Time_start"]);
  //
  uvec which_event = as<uvec>(model_data["which_event"]) - 1;
  uvec which_right = as<uvec>(model_data["which_right"]) - 1;
  uvec which_right_event = join_cols(which_event, which_right);
  uvec which_left = as<uvec>(model_data["which_left"]) - 1;
  uvec which_interval = as<uvec>(model_data["which_interval"]) - 1;
  //
  mat W0_H = as<mat>(model_data["W0_H"]);
  mat W0_h = as<mat>(model_data["W0_h"]);
  mat W0_H2 = as<mat>(model_data["W0_H2"]);
  mat W_H = as<mat>(model_data["W_H"]);
  mat W_h = as<mat>(model_data["W_h"]);
  mat W_H2 = as<mat>(model_data["W_H2"]);
  mat W_bar = as<mat>(model_data["W_bar"]);
  //
  field<mat> X_H = List2Field_mat(as<List>(model_data["X_H"]));
  field<mat> X_h = List2Field_mat(as<List>(model_data["X_h"]));
  field<mat> X_H2 = List2Field_mat(as<List>(model_data["X_H2"]));
  field<mat> Z_H = List2Field_mat(as<List>(model_data["Z_H"]));
  field<mat> Z_h = List2Field_mat(as<List>(model_data["Z_h"]));
  field<mat> Z_H2 = List2Field_mat(as<List>(model_data["Z_H2"]));
  field<mat> U_H = List2Field_mat(as<List>(model_data["U_H"]));
  field<mat> U_h = List2Field_mat(as<List>(model_data["U_h"]));
  field<mat> U_H2 = List2Field_mat(as<List>(model_data["U_H2"]));
  //
  cube vcov_prop_RE = as<cube>(vcov_prop["vcov_prop_RE"]);
  //
  mat Wlong_H = docall_cbindL(as<List>(model_data["Wlong_H"]));
  mat Wlong_h = docall_cbindL(as<List>(model_data["Wlong_h"]));
  mat Wlong_H2 = docall_cbindL(as<List>(model_data["Wlong_H2"]));
  //List Wlong_bar_ = as<List>(model_data["Wlong_bar"]);
  //field<mat> Wlong_bar = List2Field_mat(Wlong_bar_);
  field<mat> Xbase = List2Field_mat(as<List>(model_data["Xbase"]));
  // other information
  uvec idT = as<uvec>(model_data["idT"]) - 1;
  vec log_Pwk = as<vec>(model_data["log_Pwk"]);
  vec log_Pwk2 = as<vec>(model_data["log_Pwk2"]);
  uvec id_H = as<uvec>(model_data["id_H"]) - 1;
  uvec id_h = as<uvec>(model_data["id_h"]) - 1;
  uvec id_H_fast = create_fast_ind(id_H + 1);
  bool any_gammas = as<bool>(model_data["any_gammas"]);
  bool any_event = which_event.n_rows > 0;
  bool any_interval = which_interval.n_rows > 0;
  field<uvec> FunForms = List2Field_uvec(as<List>(model_info["FunForms_cpp"]), true);
  field<uvec> FunForms_ind = List2Field_uvec(as<List>(model_info["FunForms_ind"]), true);
  field<uvec> baseline = List2Field_uvec(as<List>(model_data["baseline"]), true);
  field<uvec> x_in_z = List2Field_uvec(as<List>(model_data["x_in_z"]), true);
  field<uvec> x_notin_z = List2Field_uvec(as<List>(model_data["x_notin_z"]), true);
  field<uvec> unq_idL = List2Field_uvec(as<List>(model_data["unq_idL"]), true);
  //List ind_RE_ = as<List>(model_info["ind_RE"]);
  //field<uvec> ind_RE = List2Field_uvec(ind_RE_, true);
  // initial values
  vec bs_gammas = as<vec>(initial_values["bs_gammas"]);
  vec gammas = as<vec>(initial_values["gammas"]);
  vec alphas = as<vec>(initial_values["alphas"]);
  double tau_bs_gammas = as<double>(initial_values["tau_bs_gammas"]);
  field<mat> b = List2Field_mat(as<List>(initial_values["b"]));
  mat b_mat = docall_cbindF(b);
  field<mat> mean_u(b.n_elem);
  for (uword i = 0; i < b.n_elem; ++i) mean_u.at(i) = zeros<mat>(size(b.at(i)));
  mat D = as<mat>(initial_values["D"]);
  vec sds = sqrt(D.diag());
  mat R = cov2cor(D);
  mat L = chol(R);
  field<vec> betas = List2Field_vec(as<List>(initial_values["betas"]));
  // indexes or other useful things
  uvec upper_part = trimatu_ind(size(R),  1);
  // MCMC settings
  uword n_iter = as<uword>(control["n_iter"]);
  uword n_burnin = as<uword>(control["n_burnin"]);
  bool MALA = as<bool>(control["MALA"]);
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
  double prior_D_sds_df = as<double>(priors["prior_D_sds_df"]);
  double prior_D_sds_sigma = as<double>(priors["prior_D_sds_sigma"]);
  double prior_D_L_etaLKJ = as<double>(priors["prior_D_L_etaLKJ"]);
  // store results
  uword n_b = b_mat.n_rows;
  uword n_bs_gammas = bs_gammas.n_rows;
  uword n_gammas = gammas.n_rows;
  uword n_alphas = alphas.n_rows;
  uword n_sds = sds.n_rows;
  uword n_L = vec(L(upper_part)).n_rows;
  mat res_bs_gammas(n_iter, n_bs_gammas);
  mat acceptance_bs_gammas(n_iter, n_bs_gammas, fill::zeros);
  mat res_gammas(n_iter, n_gammas);
  vec res_W_bar_gammas(n_iter);
  mat acceptance_gammas(n_iter, n_gammas, fill::zeros);
  mat res_alphas(n_iter, n_alphas);
  mat acceptance_alphas(n_iter, n_alphas, fill::zeros);
  mat res_tau_bs_gammas(n_iter, 1, fill::zeros);
  mat res_sds(n_iter, n_sds, fill::zeros);
  mat acceptance_sds(n_iter, n_sds, fill::zeros);
  mat res_L(n_iter, n_L, fill::zeros);
  mat acceptance_L(n_iter, n_L, fill::zeros);
  // scales
  vec scale_bs_gammas = create_init_scale(n_bs_gammas);
  vec scale_gammas = create_init_scale(n_gammas);
  vec scale_alphas = create_init_scale(n_alphas);
  vec scale_sds = create_init_scale(n_sds);
  vec scale_L = create_init_scale(n_L);
  vec scale_sigmas = create_init_scale(n_b);
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
  for (uword it = 0; it < n_iter; ++it) {
    double denominator_surv =
      sum(log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
               WH_gammas, Wh_gammas, WH2_gammas,
               WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
               log_Pwk, log_Pwk2, id_H_fast,
               which_event, which_right_event, which_left,
               any_interval, which_interval)) +
      logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
      logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
      logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    update_bs_gammas(bs_gammas, gammas, alphas,
                     W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                     WH_gammas, Wh_gammas, WH2_gammas,
                     WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                     log_Pwk, log_Pwk2, id_H_fast,
                     which_event, which_right_event, which_left, which_interval,
                     any_event, any_interval,
                     prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas,
                     prior_mean_gammas, prior_Tau_gammas,
                     prior_mean_alphas, prior_Tau_alphas,
                     denominator_surv, it,
                     /////
                     W0_H, W0_h, W0_H2, scale_bs_gammas, acceptance_bs_gammas,
                     res_bs_gammas);
    ////////////////////////////////////////////////////////////////////////
    double post_B_tau_bs_gammas = prior_B_tau_bs_gammas +
      0.5 * as_scalar(bs_gammas.t() * prior_Tau_bs_gammas * bs_gammas);
    tau_bs_gammas = R::rgamma(post_A_tau_bs_gammas, 1 / post_B_tau_bs_gammas);
    res_tau_bs_gammas.at(it, 0) = tau_bs_gammas;
    ////////////////////////////////////////////////////////////////////////
    if (any_gammas) {
      update_gammas(bs_gammas, gammas, alphas,
                    W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                    WH_gammas, Wh_gammas, WH2_gammas,
                    WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                    log_Pwk, log_Pwk2, id_H_fast,
                    which_event, which_right_event, which_left, which_interval,
                    any_event, any_interval,
                    prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas,
                    prior_mean_gammas, prior_Tau_gammas,
                    prior_mean_alphas, prior_Tau_alphas,
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
                  log_Pwk, log_Pwk2, id_H_fast,
                  which_event, which_right_event, which_left, which_interval,
                  any_event, any_interval,
                  prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas,
                  prior_mean_gammas, prior_Tau_gammas,
                  prior_mean_alphas, prior_Tau_alphas,
                  denominator_surv, it,
                  /////
                  Wlong_H, Wlong_h, Wlong_H2, scale_alphas,
                  acceptance_alphas, res_alphas);
    ////////////////////////////////////////////////////////////////////////
    update_D(L, sds, b_mat, upper_part,
             prior_D_sds_df, prior_D_sds_sigma, prior_D_L_etaLKJ,
             it, MALA, res_sds, res_L, scale_sds, scale_L,
             acceptance_sds, acceptance_L);
    ////////////////////////////////////////////////////////////////////////
    // update_b()...
    update_mean_u(mean_u, betas, Xbase, x_in_z, baseline, unq_idL);
    //update_Wlong(Wlong_H, Wlong_h, Wlong_H2, X_H, X_h, X_H2, Z_H, Z_h, Z_H2,
    //             U_H, U_h, U_H2, betas, b, id_H, id_h, FunForms, FunForms_ind,
    //             any_event, any_interval);
    ////////////////////////////////////////////////////////////////////////
  }
  return List::create(
    Named("mcmc") = List::create(
      Named("bs_gammas") = res_bs_gammas.rows(n_burnin, n_iter - 1),
      Named("tau_bs_gammas") = res_tau_bs_gammas.rows(n_burnin, n_iter - 1),
      Named("gammas") = res_gammas.rows(n_burnin, n_iter - 1),
      Named("W_bar_gammas") = res_W_bar_gammas.rows(n_burnin, n_iter - 1),
      Named("alphas") = res_alphas.rows(n_burnin, n_iter - 1),
      Named("sds") = res_sds.rows(n_burnin, n_iter - 1),
      Named("L") = res_L.rows(n_burnin, n_iter - 1)
    ),
    Named("acc_rate") = List::create(
      Named("bs_gammas") = acceptance_bs_gammas.rows(n_burnin, n_iter - 1),
      Named("gammas") = acceptance_gammas.rows(n_burnin, n_iter - 1),
      Named("alphas") = acceptance_alphas.rows(n_burnin, n_iter - 1),
      Named("sds") = acceptance_sds.rows(n_burnin, n_iter - 1),
      Named("L") = acceptance_L.rows(n_burnin, n_iter - 1)
    ),
    Named("mean_u") = mean_u,
    Named("Wlong_H") = Wlong_H,
    Named("Wlong_h") = Wlong_h,
    Named("Wlong_H2") = Wlong_H2
  );
}

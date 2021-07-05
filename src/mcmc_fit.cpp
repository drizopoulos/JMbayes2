#include <RcppArmadillo.h>
#include "JMbayes2_D.h"
#include "JMbayes2_Surv.h"
#include "JMbayes2_LogDens.h"
#include "JMbayes2_sigmas.h"
#include "JMbayes2_RE.h"
#include "JMbayes2_FE.h"
#include "JMbayes2_penalties.h"


// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List mcmc_cpp (List model_data, List model_info, List initial_values,
               List priors, List control) {
  // outcome vectors and design matrices
  vec Time_right = as<vec>(model_data["Time_right"]);
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
  mat W_sds = as<mat>(model_data["W_sds"]);
  mat W_std = as<mat>(model_data["W_std"]);
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
  field<mat> X = List2Field_mat(as<List>(model_data["X"]));
  field<mat> Xbar = List2Field_mat(as<List>(model_data["Xbar"]));
  field<mat> Z = List2Field_mat(as<List>(model_data["Z"]));
  field<mat> y = List2Field_mat(as<List>(model_data["y"]));
  //
  mat Wlong_H = docall_cbindL(as<List>(model_data["Wlong_H"]));
  mat Wlong_h = docall_cbindL(as<List>(model_data["Wlong_h"]));
  mat Wlong_H2 = docall_cbindL(as<List>(model_data["Wlong_H2"]));
  mat Wlong_bar = docall_cbindL(as<List>(model_data["Wlong_bar"]));
  mat Wlong_sds = docall_cbindL(as<List>(model_data["Wlong_sds"]));
  mat Wlong_std = docall_cbindL(as<List>(model_data["Wlong_std"]));
  mat X_dot = as<mat>(model_data["X_dot"]);
  // other information
  uvec idT = as<uvec>(model_data["idT"]) - 1;
  vec log_Pwk = as<vec>(model_data["log_Pwk"]);
  vec log_Pwk2 = as<vec>(model_data["log_Pwk2"]);
  uvec id_H = as<uvec>(model_data["id_H"]) - 1;
  uvec id_H_ = as<uvec>(model_data["id_H_"]) - 1;
  uvec id_h = as<uvec>(model_data["id_h"]) - 1;
  uvec id_H_fast = create_fast_ind(id_H + 1);
  uvec id_h_fast = create_fast_ind(id_h + 1);
  bool any_gammas = as<bool>(model_data["any_gammas"]);
  bool any_event = which_event.n_rows > 0;
  bool any_interval = which_interval.n_rows > 0;
  field<uvec> FunForms = List2Field_uvec(as<List>(model_info["FunForms_cpp"]), true);
  field<uvec> FunForms_ind = List2Field_uvec(as<List>(model_info["FunForms_ind"]), true);
  List Funs_FunForms = as<List>(model_info["Funs_FunForms"]);
  field<uvec> x_in_z = List2Field_uvec(as<List>(model_data["x_in_z"]), true);
  field<uvec> x_notin_z = List2Field_uvec(as<List>(model_data["x_notin_z"]), true);
  field<uvec> idL = List2Field_uvec(as<List>(model_data["idL"]), true);
  field<uvec> unq_idL = List2Field_uvec(as<List>(model_data["unq_idL"]), true);
  field<uvec> idL_lp = List2Field_uvec(as<List>(model_data["idL_lp"]), true);
  //
  field<uvec> ind_FE = List2Field_uvec(as<List>(model_data["ind_FE"]), true);
  uvec ind_FE_HC = as<uvec>(model_data["ind_FE_HC"]) - 1;
  // we do not need to subtract 1 from has_tilde_betas; logical vector
  // we will need to check if has_tilde_betas == 1
  uvec has_tilde_betas = as<uvec>(model_data["has_tilde_betas"]);
  field<uvec> ind_FE_nHC = List2Field_uvec(as<List>(model_data["ind_FE_nHC"]), true);
  uvec id_patt = as<uvec>(model_data["id_patt"]) - 1;
  field<uvec> ind_RE_patt = List2Field_uvec(as<List>(model_data["ind_RE_patt"]), true);
  field<uvec> ind_FE_patt = List2Field_uvec(as<List>(model_data["ind_FE_patt"]), true);
  //
  bool save_random_effects = as<bool>(control["save_random_effects"]);
  field<uvec> idL_lp_fast(idL_lp.n_elem);
  for (uword i = 0; i < idL_lp.n_elem; ++i) {
    idL_lp_fast.at(i) = create_fast_ind(idL_lp.at(i) + 1);
  }
  vec extra_parms = as<vec>(model_data["extra_parms"]);
  field<uvec> ind_RE = List2Field_uvec(as<List>(model_data["ind_RE"]), true);
  CharacterVector families = as<CharacterVector>(model_info["family_names"]);
  CharacterVector links = as<CharacterVector>(model_info["links"]);
  // initial values
  vec bs_gammas = as<vec>(initial_values["bs_gammas"]);
  vec tau_bs_gammas = as<vec>(initial_values["tau_bs_gammas"]);
  vec gammas = as<vec>(initial_values["gammas"]);
  vec lambda_gammas = gammas.ones();
  double tau_gammas = 1.0;
  vec nu_gammas = gammas.ones();
  double xi_gammas = 1.0;
  vec alphas = as<vec>(initial_values["alphas"]);
  vec lambda_alphas = alphas.ones();
  double tau_alphas = 1.0;
  vec nu_alphas = gammas.ones();
  double xi_alphas = 1.0;
  field<mat> b = List2Field_mat(as<List>(initial_values["b"]));
  mat b_mat = docall_cbindF(b);
  mat D = as<mat>(initial_values["D"]);
  vec sds = sqrt(D.diag());
  mat R = cov2cor(D);
  mat L = chol(R);
  field<vec> betas = List2Field_vec(as<List>(initial_values["betas"]));
  vec betas_vec = docall_rbindF(betas);
  vec sigmas = exp(as<vec>(initial_values["log_sigmas"]));
  //field<vec> sigmas = List2Field_vec(as<List>(initial_values["sigmas"]));
  uvec has_sigmas = as<uvec>(model_data["has_sigmas"]);
  // indexes or other useful things
  uvec upper_part = trimatu_ind(size(R),  1);
  uword GK_k = as<uword>(control["GK_k"]);
  umat ni_event = as<umat>(model_data["ni_event"]);
  // MCMC settings
  uword n_iter = as<uword>(control["n_iter"]);
  uword n_burnin = as<uword>(control["n_burnin"]);
  bool MALA = as<bool>(control["MALA"]);
  // priors
  field<vec> mean_bs_gammas = List2Field_vec(as<List>(priors["mean_bs_gammas"]));
  field<mat> Tau_bs_gammas = List2Field_mat(as<List>(priors["Tau_bs_gammas"]));
  vec post_A_tau_bs_gammas = as<vec>(priors["A_tau_bs_gammas"]) +
    0.5 * as<vec>(priors["rank_Tau_bs_gammas"]);
  vec B_tau_bs_gammas = as<vec>(priors["B_tau_bs_gammas"]);
  vec mean_gammas = as<vec>(priors["mean_gammas"]);
  mat Tau_gammas = as<mat>(priors["Tau_gammas"]);
  std::string penalty_gammas = as<std::string>(priors["penalty_gammas"]);
  bool shrink_gammas = penalty_gammas != "none";
  bool single_gammas = penalty_gammas == "ridge";
  double A_lambda_gammas = as<double>(priors["A_lambda_gammas"]);
  double B_lambda_gammas = as<double>(priors["B_lambda_gammas"]);
  double A_tau_gammas = as<double>(priors["A_tau_gammas"]);
  double B_tau_gammas = as<double>(priors["B_tau_gammas"]);
  double A_nu_gammas = as<double>(priors["A_nu_gammas"]);
  double B_nu_gammas = as<double>(priors["B_nu_gammas"]);
  double A_xi_gammas = as<double>(priors["A_xi_gammas"]);
  double B_xi_gammas = as<double>(priors["B_xi_gammas"]);
  vec mean_alphas = as<vec>(priors["mean_alphas"]);
  mat Tau_alphas = as<mat>(priors["Tau_alphas"]);
  std::string penalty_alphas = as<std::string>(priors["penalty_alphas"]);
  bool shrink_alphas = penalty_alphas != "none";
  bool single_alphas = penalty_alphas == "ridge";
  double A_lambda_alphas = as<double>(priors["A_lambda_alphas"]);
  double B_lambda_alphas = as<double>(priors["B_lambda_alphas"]);
  double A_tau_alphas = as<double>(priors["A_tau_alphas"]);
  double B_tau_alphas = as<double>(priors["B_tau_alphas"]);
  double A_nu_alphas = as<double>(priors["A_nu_alphas"]);
  double B_nu_alphas = as<double>(priors["B_nu_alphas"]);
  double A_xi_alphas = as<double>(priors["A_xi_alphas"]);
  double B_xi_alphas = as<double>(priors["B_xi_alphas"]);
  bool gamma_prior_D_sds = as<bool>(priors["gamma_prior_D_sds"]);
  double D_sds_df = as<double>(priors["D_sds_df"]);
  vec D_sds_sigma = as<vec>(priors["D_sds_sigma"]);
  double D_sds_shape = as<double>(priors["D_sds_shape"]);
  vec D_sds_mean = as<vec>(priors["D_sds_mean"]);
  double D_L_etaLKJ = as<double>(priors["D_L_etaLKJ"]);
  bool gamma_prior_sigmas = as<bool>(priors["gamma_prior_sigmas"]);
  double sigmas_df = as<double>(priors["sigmas_df"]);
  vec sigmas_sigmas = as<vec>(priors["sigmas_sigmas"]);
  double sigmas_shape = as<double>(priors["sigmas_shape"]);
  vec sigmas_mean = as<vec>(priors["sigmas_mean"]);
  vec mean_betas_HC = as<vec>(priors["mean_betas_HC"]);
  mat Tau_betas_HC = as<mat>(priors["Tau_betas_HC"]);
  vec Tau_mean_betas_HC = Tau_betas_HC * mean_betas_HC;
  field<vec> mean_betas_nHC = List2Field_vec(as<List>(priors["mean_betas_nHC"]));
  field<mat> Tau_betas_nHC = List2Field_mat(as<List>(priors["Tau_betas_nHC"]));
  // store results
  uword n_b = b_mat.n_rows;
  uword n_bs_gammas = bs_gammas.n_rows;
  uword n_strata = tau_bs_gammas.n_rows;
  uword n_per_stratum = n_bs_gammas / n_strata;
  uword n_gammas = gammas.n_rows;
  uword n_alphas = alphas.n_rows;
  uword n_sds = sds.n_rows;
  uword n_L = vec(L(upper_part)).n_rows;
  uword n_sigmas = sigmas.n_rows;
  // uword n_sigmas = n_field(sigmas);
  uword n_betas = betas_vec.n_rows;
  mat res_bs_gammas(n_iter, n_bs_gammas, fill::zeros);
  mat acceptance_bs_gammas(n_iter, n_bs_gammas, fill::zeros);
  mat res_gammas(n_iter, n_gammas, fill::zeros);
  mat res_W_std_gammas(n_iter, 1, fill::zeros);
  mat acceptance_gammas(n_iter, n_gammas, fill::zeros);
  mat res_alphas(n_iter, n_alphas, fill::zeros);
  mat res_Wlong_std_alphas(n_iter, 1, fill::zeros);
  mat acceptance_alphas(n_iter, n_alphas, fill::zeros);
  mat res_tau_bs_gammas(n_iter, n_strata, fill::zeros);
  mat res_sds(n_iter, n_sds, fill::zeros);
  mat acceptance_sds(n_iter, n_sds, fill::zeros);
  mat res_L(n_iter, n_L, fill::zeros);
  mat acceptance_L(n_iter, n_L, fill::zeros);
  cube res_b(n_b, b_mat.n_cols, 1, fill::zeros);
  cube res_b_last(n_b, b_mat.n_cols, 1, fill::zeros);
  if (save_random_effects) {
    res_b.set_size(n_b, b_mat.n_cols, n_iter);
  }
  mat acceptance_b(n_b, b_mat.n_cols, fill::zeros);
  mat cumsum_b(n_b, b_mat.n_cols);
  cube outprod_b(b_mat.n_cols, b_mat.n_cols, b_mat.n_rows);
  cube var_b(b_mat.n_cols, b_mat.n_cols, b_mat.n_rows);
  mat res_sigmas(n_iter, n_sigmas, fill::zeros);
  mat acceptance_sigmas(n_iter, n_sigmas, fill::zeros);
  mat res_betas(n_iter, n_betas, fill::zeros);
  field<vec> acceptance_betas = create_storage(x_notin_z);
  mat res_logLik(n_iter, n_b, fill::zeros);
  // scales
  vec scale_bs_gammas = create_init_scale(n_bs_gammas);
  vec scale_gammas = create_init_scale(n_gammas);
  vec scale_alphas = create_init_scale(n_alphas);
  vec scale_sds = create_init_scale(n_sds);
  vec scale_L = create_init_scale(n_L);
  mat scale_b = mat(n_b,  b_mat.n_cols, fill::ones) * 0.1;
  vec scale_sigmas = create_init_scale(n_sigmas);
  field<vec> scale_betas = create_init_scaleF(x_notin_z);
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
    WH2_gammas = W_H2 * gammas;
  }
  vec WlongH_alphas = Wlong_H * alphas;
  vec Wlongh_alphas(W0_h.n_rows);
  vec WlongH2_alphas(W0_H2.n_rows);
  if (any_event) {
    Wlongh_alphas = Wlong_h * alphas;
  }
  if (any_interval) {
    WlongH2_alphas = Wlong_H2 * alphas;
  }
  //
  uvec which_term_h = as<uvec>(model_data["which_term_h"]) - 1; //!! new //?? later split these elements in the different sections
  uvec which_term_H = as<uvec>(model_data["which_term_H"]) - 1; //!! new
  bool recurrent = as<bool>(model_data["recurrent"]); //!! new
  vec alphaF = as<vec>(initial_values["alphaF"]); //!! new
  vec frailty = as<vec>(initial_values["frailty"]); //!! new
  alphaF.zeros(); //?? delete later
  frailty.zeros(); //?? delete later 
  mat res_alphaF(n_iter, 1, fill::zeros); //!! new
  mat acceptance_alphaF(n_iter, 1, fill::zeros); //!! new
  mat res_frailty(n_iter, frailty.n_elem, fill::zeros); //!! new
  mat acceptance_frailty(n_iter, frailty.n_elem, fill::zeros); //!! new
  vec alphaF_H(WH_gammas.n_rows, fill::ones); //!! new
  vec alphaF_h(Wh_gammas.n_rows, fill::ones); //!! new
  alphaF_H.rows(which_term_H).fill(alphaF.at(1)); //!! new
  alphaF_h.rows(which_term_h).fill(alphaF.at(1)); //!! new
  vec frailty_H(WH_gammas.n_rows, fill::zeros); //!! new
  vec frailty_h(Wh_gammas.n_rows, fill::zeros); //!! new
  frailty_h = frailty.rows(id_h); //!! new
  frailty_H = frailty.rows(id_H_); //!! new
  //
  vec logLik_surv =
    log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
             WH_gammas, Wh_gammas, WH2_gammas,
             WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
             log_Pwk, log_Pwk2, id_H_fast, id_h_fast,
             which_event, which_right_event, which_left,
             any_interval, which_interval,
             recurrent, frailty_H, frailty_h, alphaF_H, alphaF_h); //!! new
  double denominator_surv =
    sum(logLik_surv) +
    logPrior_surv(bs_gammas, gammas, alphas, mean_bs_gammas,
                  Tau_bs_gammas, tau_bs_gammas,
                  mean_gammas, Tau_gammas, lambda_gammas,
                  tau_gammas, shrink_gammas,
                  mean_alphas, Tau_alphas, lambda_alphas,
                  tau_alphas, shrink_alphas);
  //
  vec logLik_re = log_re(b_mat, L, sds);
  //
  field<vec> eta = linpred_mixed(X, betas, Z, b, idL);
  vec logLik_long = log_long(y, eta, sigmas, extra_parms, families, links,
                             idL_lp_fast, unq_idL, n_b);
  //
  for (uword it = 0; it < n_iter; ++it) {
    update_bs_gammas(bs_gammas, gammas, alphas,
                     W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                     WH_gammas, Wh_gammas, WH2_gammas,
                     WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                     log_Pwk, log_Pwk2, id_H_fast, id_h_fast,
                     which_event, which_right_event, which_left, which_interval,
                     any_event, any_interval,
                     mean_bs_gammas, Tau_bs_gammas, tau_bs_gammas,
                     mean_gammas, Tau_gammas, lambda_gammas,
                     tau_gammas, shrink_gammas,
                     mean_alphas, Tau_alphas, lambda_alphas,
                     tau_alphas, shrink_alphas,
                     logLik_surv, denominator_surv, it,
                     /////
                     W0_H, W0_h, W0_H2, scale_bs_gammas, acceptance_bs_gammas,
                     res_bs_gammas,
                     recurrent, frailty_H, frailty_h, alphaF_H, alphaF_h //!! new
                    );

    ////////////////////////////////////////////////////////////////////////

    for (uword j = 0; j < n_strata; ++j) {
      vec bs_gammas_j =
        bs_gammas.rows(j * n_per_stratum, (j + 1) * n_per_stratum - 1);
      double quad = as_scalar(bs_gammas_j.t() * Tau_bs_gammas.at(j) *
                              bs_gammas_j);
      double post_B_tau = B_tau_bs_gammas.at(j) + 0.5 * quad;
      tau_bs_gammas.at(j) = R::rgamma(post_A_tau_bs_gammas.at(j), 1 / post_B_tau);
      res_tau_bs_gammas.at(it, j) = tau_bs_gammas.at(j);
    }

    ////////////////////////////////////////////////////////////////////////

    if (any_gammas) {
      update_gammas(bs_gammas, gammas, alphas,
                    W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                    WH_gammas, Wh_gammas, WH2_gammas,
                    WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                    log_Pwk, log_Pwk2, id_H_fast, id_h_fast,
                    which_event, which_right_event, which_left, which_interval,
                    any_event, any_interval,
                    mean_bs_gammas, Tau_bs_gammas, tau_bs_gammas,
                    mean_gammas, Tau_gammas, lambda_gammas,
                    tau_gammas, shrink_gammas,
                    mean_alphas, Tau_alphas, lambda_alphas,
                    tau_alphas, shrink_alphas,
                    logLik_surv, denominator_surv, it,
                    /////
                    W_H, W_h, W_H2, scale_gammas, acceptance_gammas,
                    res_gammas,
                    //
                    recurrent, frailty_H, frailty_h, alphaF_H, alphaF_h //!! new
                    );
      res_W_std_gammas.at(it) = as_scalar(W_std * gammas);

      if (shrink_gammas) {
        update_penalties (
            gammas, lambda_gammas, tau_gammas, nu_gammas, xi_gammas,
            single_gammas, A_lambda_gammas, B_lambda_gammas,
            A_tau_gammas, B_tau_gammas, A_nu_gammas, B_nu_gammas,
            A_xi_gammas, B_xi_gammas);
      }
    }

    ////////////////////////////////////////////////////////////////////////

    update_alphas(bs_gammas, gammas, alphas,
                  W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                  WH_gammas, Wh_gammas, WH2_gammas,
                  WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                  log_Pwk, log_Pwk2, id_H_fast, id_h_fast,
                  which_event, which_right_event, which_left, which_interval,
                  any_event, any_interval,
                  mean_bs_gammas, Tau_bs_gammas, tau_bs_gammas,
                  mean_gammas, Tau_gammas, lambda_gammas,
                  tau_gammas, shrink_gammas,
                  mean_alphas, Tau_alphas, lambda_alphas,
                  tau_alphas, shrink_alphas,
                  logLik_surv, denominator_surv, it,
                  /////
                  Wlong_H, Wlong_h, Wlong_H2, scale_alphas,
                  acceptance_alphas, res_alphas,
                  //
                  recurrent, frailty_H, frailty_h, alphaF_H, alphaF_h //!! new
                  );

    res_Wlong_std_alphas.at(it) = as_scalar(Wlong_std * alphas);

    if (shrink_alphas) {
      update_penalties (
          alphas, lambda_alphas, tau_alphas, nu_alphas, xi_alphas,
          single_alphas, A_lambda_alphas, B_lambda_alphas,
          A_tau_alphas, B_tau_alphas, A_nu_alphas, B_nu_alphas,
          A_xi_alphas, B_xi_alphas);
    }

    ////////////////////////////////////////////////////////////////////////

    if (it > 99) {

    update_D(L, sds, b_mat, upper_part,
             D_sds_df, D_sds_sigma, D_sds_shape, D_sds_mean, gamma_prior_D_sds,
             D_L_etaLKJ,
             it, MALA, logLik_re, res_sds, res_L, scale_sds, scale_L,
             acceptance_sds, acceptance_L);

    ////////////////////////////////////////////////////////////////////////

    update_b(b, b_mat, eta, logLik_long, logLik_surv, logLik_re,
             Wlong_H, Wlong_h, Wlong_H2, WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
             scale_b, ind_RE,
             X_H, X_h, X_H2, Z_H, Z_h, Z_H2, U_H, U_h, U_H2,
             Wlong_bar, Wlong_sds, betas, alphas, id_H_, id_h,
             FunForms, FunForms_ind, Funs_FunForms, X, Z, idL, y, sigmas,
             extra_parms, families, links, idL_lp_fast, unq_idL,
             W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas, WH_gammas,
             Wh_gammas, WH2_gammas, log_Pwk, log_Pwk2,
             id_H_fast, id_h_fast, which_event, which_right_event, which_left,
             which_interval, any_event, any_interval, ni_event,
             L, sds, it, acceptance_b, res_b, res_b_last, save_random_effects,
             n_burnin, n_iter, GK_k, cumsum_b, outprod_b,
             //
             recurrent, frailty_H, frailty_h, alphaF_H, alphaF_h //!! new
             );

    ////////////////////////////////////////////////////////////////////

    update_sigmas(sigmas, has_sigmas, y, eta, extra_parms, families, links,
                  idL_lp_fast, gamma_prior_sigmas, sigmas_df, sigmas_sigmas,
                  sigmas_shape, sigmas_mean, it, res_sigmas, scale_sigmas,
                  acceptance_sigmas);

    logLik_long = log_long(y, eta, sigmas, extra_parms, families, links,
                           idL_lp_fast, unq_idL, n_b);

    ////////////////////////////////////////////////////////////////////

    update_betas(betas, res_betas, acceptance_betas, scale_betas, eta,
                 logLik_long, logLik_surv, Wlong_H, Wlong_h, Wlong_H2,
                 WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                 Tau_mean_betas_HC, Tau_betas_HC, b_mat, L, sds, X_dot,
                 ind_FE, ind_FE_HC, id_patt, ind_RE_patt, ind_FE_patt,
                 it, has_tilde_betas, X, Z, b, idL, y, sigmas,
                 extra_parms, families, links, idL_lp_fast, mean_betas_nHC,
                 Tau_betas_nHC, x_notin_z,
                 X_H, X_h, X_H2, Z_H, Z_h, Z_H2, U_H, U_h, U_H2,
                 Wlong_bar, Wlong_sds, id_H_, id_h, FunForms, FunForms_ind,
                 Funs_FunForms, alphas, any_event, any_interval,
                 W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                 WH_gammas, Wh_gammas, WH2_gammas,
                 log_Pwk, log_Pwk2, id_H_fast, id_h_fast, which_event,
                 which_right_event, which_left, which_interval, unq_idL, n_burnin,
                 //
                 recurrent, frailty_H, frailty_h, alphaF_H, alphaF_h //!! new
                 );

    // update intercepts
    for (uword j = 0; j < y.n_elem; ++j) {
      res_betas.at(it, ind_FE.at(j).at(0)) -= as_scalar(Xbar.at(j) * betas.at(j));
    }

    denominator_surv =
      sum(logLik_surv) +
      logPrior_surv(bs_gammas, gammas, alphas, mean_bs_gammas,
                    Tau_bs_gammas, tau_bs_gammas,
                    mean_gammas, Tau_gammas, lambda_gammas,
                    tau_gammas, shrink_gammas,
                    mean_alphas, Tau_alphas, lambda_alphas,
                    tau_alphas, shrink_alphas);
    }
    ////////////////////////////////////////////////////////////////////

    res_logLik.row(it) = trans(logLik_long + logLik_surv + logLik_re);
  }
  if (save_random_effects) {
    res_b = res_b.slices(n_burnin, n_iter - 1);
  } else {
    res_b.slice(0) = cumsum_b / (n_iter - n_burnin);
  }
  acceptance_b = acceptance_b / (n_iter - n_burnin);
  for (uword j = 0; j < acceptance_betas.n_elem; ++j) {
    acceptance_betas.at(j) = acceptance_betas.at(j) / (n_iter - n_burnin);
  }
  if (any_gammas) {
    res_gammas.each_row() /= W_sds;
  }
  res_alphas.each_row() /= Wlong_sds;
  return List::create(
    Named("mcmc") = List::create(
      Named("bs_gammas") = res_bs_gammas.rows(n_burnin, n_iter - 1),
      Named("tau_bs_gammas") = res_tau_bs_gammas.rows(n_burnin, n_iter - 1),
      Named("gammas") = res_gammas.rows(n_burnin, n_iter - 1),
      Named("W_std_gammas") = res_W_std_gammas.rows(n_burnin, n_iter - 1),
      Named("alphas") = res_alphas.rows(n_burnin, n_iter - 1),
      Named("Wlong_std_alphas") = res_Wlong_std_alphas.rows(n_burnin, n_iter - 1),
      Named("sds") = res_sds.rows(n_burnin, n_iter - 1),
      Named("L") = res_L.rows(n_burnin, n_iter - 1),
      Named("b") = res_b,
      Named("b_last") = res_b_last,
      Named("cumsum_b") = cumsum_b,
      Named("outprod_b") = outprod_b,
      Named("sigmas") = res_sigmas.rows(n_burnin, n_iter - 1),
      Named("betas") = res_betas.rows(n_burnin, n_iter - 1),
      Named("alphaF") = res_alphaF.rows(n_burnin, n_iter - 1),
      Named("frailty") = res_frailty.rows(n_burnin, n_iter - 1)
    ),
    Named("acc_rate") = List::create(
      Named("bs_gammas") = mean(acceptance_bs_gammas.rows(n_burnin, n_iter - 1)),
      Named("gammas") = mean(acceptance_gammas.rows(n_burnin, n_iter - 1)),
      Named("alphas") = mean(acceptance_alphas.rows(n_burnin, n_iter - 1)),
      Named("sds") = mean(acceptance_sds.rows(n_burnin, n_iter - 1)),
      Named("L") = mean(acceptance_L.rows(n_burnin, n_iter - 1)),
      Named("b") = acceptance_b,
      Named("sigmas") = mean(acceptance_sigmas.rows(n_burnin, n_iter - 1)),
      Named("betas") = acceptance_betas,
      Named("alphaF") = mean(acceptance_alphaF.rows(n_burnin, n_iter - 1)),
      Named("frailty") = acceptance_frailty
    ),
    Named("logLik") = res_logLik.rows(n_burnin, n_iter - 1)
  );
}

// [[Rcpp::export]]
arma::vec logLik_jm (List thetas, List model_data, List model_info,
                     List control) {
  field<vec> betas = List2Field_vec(as<List>(thetas["betas"]));
  mat b_mat = as<mat>(thetas["b"]);
  field<mat> b =
    mat2field(b_mat, List2Field_uvec(as<List>(model_data["ind_RE"]), true));
  vec sigmas = as<vec>(thetas["sigmas"]);
  vec bs_gammas = as<vec>(thetas["bs_gammas"]);
  vec gammas = as<vec>(thetas["gammas"]);
  vec alphas = as<vec>(thetas["alphas"]);
  vec tau_bs_gammas = as<vec>(thetas["tau_bs_gammas"]);
  mat D = as<mat>(thetas["D"]);
  vec sds = sqrt(D.diag());
  mat R = cov2cor(D);
  mat L = chol(R);
  /////////////
  field<mat> y = List2Field_mat(as<List>(model_data["y"]));
  field<mat> X = List2Field_mat(as<List>(model_data["X"]));
  field<mat> Xbar = List2Field_mat(as<List>(model_data["Xbar"]));
  field<mat> Z = List2Field_mat(as<List>(model_data["Z"]));
  vec extra_parms = as<vec>(model_data["extra_parms"]);
  CharacterVector families = as<CharacterVector>(model_info["family_names"]);
  CharacterVector links = as<CharacterVector>(model_info["links"]);
  field<uvec> idL = List2Field_uvec(as<List>(model_data["idL"]), true);
  field<uvec> idL_lp = List2Field_uvec(as<List>(model_data["idL_lp"]), true);
  field<uvec> idL_lp_fast(idL_lp.n_elem);
  for (uword i = 0; i < idL_lp.n_elem; ++i) {
    idL_lp_fast.at(i) = create_fast_ind(idL_lp.at(i) + 1);
  }
  field<uvec> unq_idL = List2Field_uvec(as<List>(model_data["unq_idL"]), true);
  /////////////
  mat W0_H = as<mat>(model_data["W0_H"]);
  mat W0_h = as<mat>(model_data["W0_h"]);
  mat W0_H2 = as<mat>(model_data["W0_H2"]);
  mat W_H = as<mat>(model_data["W_H"]);
  mat W_h = as<mat>(model_data["W_h"]);
  mat W_H2 = as<mat>(model_data["W_H2"]);
  field<mat> X_H = List2Field_mat(as<List>(model_data["X_H"]));
  field<mat> X_h = List2Field_mat(as<List>(model_data["X_h"]));
  field<mat> X_H2 = List2Field_mat(as<List>(model_data["X_H2"]));
  field<mat> Z_H = List2Field_mat(as<List>(model_data["Z_H"]));
  field<mat> Z_h = List2Field_mat(as<List>(model_data["Z_h"]));
  field<mat> Z_H2 = List2Field_mat(as<List>(model_data["Z_H2"]));
  field<mat> U_H = List2Field_mat(as<List>(model_data["U_H"]));
  field<mat> U_h = List2Field_mat(as<List>(model_data["U_h"]));
  field<mat> U_H2 = List2Field_mat(as<List>(model_data["U_H2"]));
  mat Wlong_bar = docall_cbindL(as<List>(model_data["Wlong_bar"]));
  mat Wlong_sds = docall_cbindL(as<List>(model_data["Wlong_sds"]));
  mat W_sds = as<mat>(model_data["W_sds"]);
  uvec which_event = as<uvec>(model_data["which_event"]) - 1;
  uvec which_right = as<uvec>(model_data["which_right"]) - 1;
  uvec which_right_event = join_cols(which_event, which_right);
  uvec which_left = as<uvec>(model_data["which_left"]) - 1;
  uvec which_interval = as<uvec>(model_data["which_interval"]) - 1;
  bool any_gammas = as<bool>(model_data["any_gammas"]);
  bool any_event = which_event.n_rows > 0;
  bool any_interval = which_interval.n_rows > 0;
  field<uvec> FunForms = List2Field_uvec(as<List>(model_info["FunForms_cpp"]), true);
  field<uvec> FunForms_ind = List2Field_uvec(as<List>(model_info["FunForms_ind"]), true);
  List Funs_FunForms = as<List>(model_info["Funs_FunForms"]);
  uvec id_H_ = as<uvec>(model_data["id_H_"]) - 1;
  uvec id_h = as<uvec>(model_data["id_h"]) - 1;
  vec log_Pwk = as<vec>(model_data["log_Pwk"]);
  vec log_Pwk2 = as<vec>(model_data["log_Pwk2"]);
  uvec id_H = as<uvec>(model_data["id_H"]) - 1;
  uvec id_H_fast = create_fast_ind(id_H + 1);
  uvec id_h_fast = create_fast_ind(id_h + 1);
  //
  uvec which_term_h = as<uvec>(model_data["which_term_h"]) - 1; //!! new
  uvec which_term_H = as<uvec>(model_data["which_term_H"]) - 1; //!! new
  bool recurrent = as<bool>(model_data["recurrent"]); //!! new
  vec alphaF = as<vec>(thetas["alphaF"]); //!! new
  vec frailty = as<vec>(thetas["frailty"]); //!! new
  /////////////
  vec out =
    logLik_jm_stripped(
      betas, b, sigmas, bs_gammas, gammas, alphas, tau_bs_gammas, L, sds,
      ///
      y, X, Xbar, Z, extra_parms, families, links, idL, idL_lp_fast, unq_idL,
      ///
      W0_H, W0_h, W0_H2, W_H, W_h, W_H2, X_H, X_h, X_H2, Z_H, Z_h, Z_H2,
      U_H, U_h, U_H2, Wlong_bar, Wlong_sds, W_sds, any_event, any_interval, any_gammas,
      FunForms, FunForms_ind, Funs_FunForms, id_H_, id_h, log_Pwk, log_Pwk2,
      id_H_fast, id_h_fast, which_event, which_right_event, which_left,
      which_interval,
      recurrent, alphaF, frailty, which_term_H, which_term_h //!! new
    );
  return out;
}

// [[Rcpp::export]]
arma::mat mlogLik_jm (List res_thetas, arma::mat mean_b_mat, arma::cube post_vars,
                      List model_data, List model_info, List control) {
  field<mat> betas = List2Field_mat(as<List>(res_thetas["betas"]));
  for (uword j = 0; j < betas.n_elem; ++j) betas.at(j) = trans(betas.at(j));
  mat sigmas = trans(as<mat>(res_thetas["sigmas"]));
  mat bs_gammas = trans(as<mat>(res_thetas["bs_gammas"]));
  mat gammas = trans(as<mat>(res_thetas["gammas"]));
  mat alphas = trans(as<mat>(res_thetas["alphas"]));
  mat tau_bs_gammas = trans(as<mat>(res_thetas["tau_bs_gammas"]));
  cube D = as<cube>(res_thetas["D"]);
  uword K = D.n_slices;
  uword q = D.slice(0).n_rows;
  mat sds(q, K);
  cube L(q, q, K);
  for (uword i = 0; i < K; ++i) {
    mat D_i = D.slice(i);
    sds.col(i) = sqrt(D_i.diag());
    mat R = cov2cor(D_i);
    L.slice(i) = chol(R);
  }
  /////////////
  uword n = mean_b_mat.n_rows;
  field<mat> mean_b =
    mat2field(mean_b_mat, List2Field_uvec(as<List>(model_data["ind_RE"]), true));
  vec log_det_post_vars(n);
  for (uword i = 0; i < n; ++i) {
    log_det_post_vars.at(i) = log_det_sympd(post_vars.slice(i));
  }
  /////////////
  field<mat> y = List2Field_mat(as<List>(model_data["y"]));
  field<mat> X = List2Field_mat(as<List>(model_data["X"]));
  field<mat> Xbar = List2Field_mat(as<List>(model_data["Xbar"]));
  field<mat> Z = List2Field_mat(as<List>(model_data["Z"]));
  vec extra_parms = as<vec>(model_data["extra_parms"]);
  CharacterVector families = as<CharacterVector>(model_info["family_names"]);
  CharacterVector links = as<CharacterVector>(model_info["links"]);
  field<uvec> idL = List2Field_uvec(as<List>(model_data["idL"]), true);
  field<uvec> idL_lp = List2Field_uvec(as<List>(model_data["idL_lp"]), true);
  field<uvec> idL_lp_fast(idL_lp.n_elem);
  for (uword i = 0; i < idL_lp.n_elem; ++i) {
    idL_lp_fast.at(i) = create_fast_ind(idL_lp.at(i) + 1);
  }
  field<uvec> unq_idL = List2Field_uvec(as<List>(model_data["unq_idL"]), true);
  field<uvec> ind_FE = List2Field_uvec(as<List>(model_data["ind_FE"]), true);
  /////////////
  mat W0_H = as<mat>(model_data["W0_H"]);
  mat W0_h = as<mat>(model_data["W0_h"]);
  mat W0_H2 = as<mat>(model_data["W0_H2"]);
  mat W_H = as<mat>(model_data["W_H"]);
  mat W_h = as<mat>(model_data["W_h"]);
  mat W_H2 = as<mat>(model_data["W_H2"]);
  field<mat> X_H = List2Field_mat(as<List>(model_data["X_H"]));
  field<mat> X_h = List2Field_mat(as<List>(model_data["X_h"]));
  field<mat> X_H2 = List2Field_mat(as<List>(model_data["X_H2"]));
  field<mat> Z_H = List2Field_mat(as<List>(model_data["Z_H"]));
  field<mat> Z_h = List2Field_mat(as<List>(model_data["Z_h"]));
  field<mat> Z_H2 = List2Field_mat(as<List>(model_data["Z_H2"]));
  field<mat> U_H = List2Field_mat(as<List>(model_data["U_H"]));
  field<mat> U_h = List2Field_mat(as<List>(model_data["U_h"]));
  field<mat> U_H2 = List2Field_mat(as<List>(model_data["U_H2"]));
  mat Wlong_bar = docall_cbindL(as<List>(model_data["Wlong_bar"]));
  mat Wlong_sds = docall_cbindL(as<List>(model_data["Wlong_sds"]));
  mat W_sds = as<mat>(model_data["W_sds"]);
  uvec which_event = as<uvec>(model_data["which_event"]) - 1;
  uvec which_right = as<uvec>(model_data["which_right"]) - 1;
  uvec which_right_event = join_cols(which_event, which_right);
  uvec which_left = as<uvec>(model_data["which_left"]) - 1;
  uvec which_interval = as<uvec>(model_data["which_interval"]) - 1;
  bool any_gammas = as<bool>(model_data["any_gammas"]);
  bool any_event = which_event.n_rows > 0;
  bool any_interval = which_interval.n_rows > 0;
  field<uvec> FunForms = List2Field_uvec(as<List>(model_info["FunForms_cpp"]), true);
  field<uvec> FunForms_ind = List2Field_uvec(as<List>(model_info["FunForms_ind"]), true);
  List Funs_FunForms = as<List>(model_info["Funs_FunForms"]);
  uvec id_H_ = as<uvec>(model_data["id_H_"]) - 1;
  uvec id_h = as<uvec>(model_data["id_h"]) - 1;
  vec log_Pwk = as<vec>(model_data["log_Pwk"]);
  vec log_Pwk2 = as<vec>(model_data["log_Pwk2"]);
  uvec id_H = as<uvec>(model_data["id_H"]) - 1;
  uvec id_H_fast = create_fast_ind(id_H + 1);
  uvec id_h_fast = create_fast_ind(id_h + 1);
  //
  uvec which_term_h = as<uvec>(model_data["which_term_h"]) - 1; //!! new
  uvec which_term_H = as<uvec>(model_data["which_term_H"]) - 1; //!! new
  bool recurrent = as<bool>(model_data["recurrent"]); //!! new
  mat alphaF = trans(as<mat>(res_thetas["alphaF"])); //!! new
  mat frailty = trans(as<mat>(res_thetas["frailty"])); //!! new
  /////////////
  mat out(n, K);
  field<vec> betas_i(betas.n_elem);
  for (uword i = 0; i < K; ++i) {
    for (uword j = 0; j < betas.n_elem; ++j) betas_i.at(j) = betas.at(j).col(i);
    vec oo = logLik_jm_stripped(
      betas_i, mean_b, sigmas.col(i), bs_gammas.col(i), gammas.col(i),
      alphas.col(i), tau_bs_gammas.col(i), L.slice(i), sds.col(i),
      ///
      y, X, Xbar, Z, extra_parms, families, links, idL, idL_lp_fast, unq_idL,
      ///
      W0_H, W0_h, W0_H2, W_H, W_h, W_H2, X_H, X_h, X_H2, Z_H, Z_h, Z_H2,
      U_H, U_h, U_H2, Wlong_bar, Wlong_sds, W_sds, any_event, any_interval, any_gammas,
      FunForms, FunForms_ind, Funs_FunForms, id_H_, id_h, log_Pwk, log_Pwk2,
      id_H_fast, id_h_fast, which_event, which_right_event, which_left,
      which_interval,
      recurrent, alphaF.col(i), frailty.col(i), which_term_H, which_term_h //!! new
    );
    oo += 0.5 * ((double)mean_b_mat.n_cols * log2pi + log_det_post_vars);
    out.col(i) = oo;
  }
  out = out.t();
  return out;
}

// [[Rcpp::export]]
arma::cube simulate_REs (List Data, List MCMC, List control) {
  //////////////////////////////
  // Event Process Components //
  //////////////////////////////
  uvec which_event = as<uvec>(Data["which_event"]) - 1;
  uvec which_right = as<uvec>(Data["which_right"]) - 1;
  uvec which_right_event = join_cols(which_event, which_right);
  uvec which_left = as<uvec>(Data["which_left"]) - 1;
  uvec which_interval = as<uvec>(Data["which_interval"]) - 1;
  bool any_event = which_event.n_rows > 0;
  bool any_interval = which_interval.n_rows > 0;
  umat ni_event = as<umat>(Data["ni_event"]);
  uvec idT = as<uvec>(Data["idT"]) - 1;
  vec log_Pwk = as<vec>(Data["log_Pwk"]);
  vec log_Pwk2 = as<vec>(Data["log_Pwk2"]);
  uvec id_H = as<uvec>(Data["id_H"]) - 1;
  uvec id_H_ = as<uvec>(Data["id_H_"]) - 1;
  uvec id_h = as<uvec>(Data["id_h"]) - 1;
  uvec indFast_H = create_fast_ind(id_H + 1);
  uvec indFast_h;
  if (id_h.n_rows> 1) indFast_h = create_fast_ind(id_h + 1); else indFast_h = id_h;
  uword GK_k = as<uword>(control["GK_k"]);
  bool any_gammas = as<bool>(Data["any_gammas"]);
  field<uvec> FunForms = List2Field_uvec(as<List>(Data["FunForms_cpp"]), true);
  field<uvec> FunForms_ind = List2Field_uvec(as<List>(Data["FunForms_ind"]), true);
  List Funs_FunForms = as<List>(Data["Funs_FunForms"]);
  //
  field<uvec> ind_RE = List2Field_uvec(as<List>(Data["ind_RE"]), true);
  mat W0_H = as<mat>(Data["W0_H"]);
  mat W0_h = as<mat>(Data["W0_h"]);
  mat W0_H2 = as<mat>(Data["W0_H2"]);
  mat W_H = as<mat>(Data["W_H"]);
  mat W_h = as<mat>(Data["W_h"]);
  mat W_H2 = as<mat>(Data["W_H2"]);
  field<mat> X_H = List2Field_mat(as<List>(Data["X_H"]));
  field<mat> X_h = List2Field_mat(as<List>(Data["X_h"]));
  field<mat> X_H2 = List2Field_mat(as<List>(Data["X_H2"]));
  field<mat> Z_H = List2Field_mat(as<List>(Data["Z_H"]));
  field<mat> Z_h = List2Field_mat(as<List>(Data["Z_h"]));
  field<mat> Z_H2 = List2Field_mat(as<List>(Data["Z_H2"]));
  field<mat> U_H = List2Field_mat(as<List>(Data["U_H"]));
  field<mat> U_h = List2Field_mat(as<List>(Data["U_h"]));
  field<mat> U_H2 = List2Field_mat(as<List>(Data["U_H2"]));
  mat Wlong_bar = docall_cbindL(as<List>(Data["Wlong_bar"]));
  mat Wlong_sds = docall_cbindL(as<List>(Data["Wlong_sds"]));
  /////////////////////////////////////
  // Longitudinal Process Components //
  /////////////////////////////////////
  field<mat> y = List2Field_mat(as<List>(Data["y"]));
  field<mat> X = List2Field_mat(as<List>(Data["X"]));
  field<mat> Z = List2Field_mat(as<List>(Data["Z"]));
  //
  vec extra_parms = as<vec>(Data["extra_parms"]);
  CharacterVector families = as<CharacterVector>(Data["family_names"]);
  CharacterVector links = as<CharacterVector>(Data["links"]);
  field<uvec> idL = List2Field_uvec(as<List>(Data["idL"]), true);
  field<uvec> unq_idL = List2Field_uvec(as<List>(Data["unq_idL"]), true);
  field<uvec> idL_lp = List2Field_uvec(as<List>(Data["idL_lp"]), true);
  field<uvec> ids(idL_lp.n_elem);
  for (uword i = 0; i < idL_lp.n_elem; ++i) {
    ids.at(i) = create_fast_ind(idL_lp.at(i) + 1);
  }
  ////////////////////////////
  // MCMC Sample Parameters //
  ///////////////////////////
  field<mat> b = List2Field_mat(as<List>(MCMC["b"]));
  mat b_mat = docall_cbindF(b);
  mat bs_gammas = trans(as<mat>(MCMC["bs_gammas"]));
  mat gammas = trans(as<mat>(MCMC["gammas"]));
  mat alphas = trans(as<mat>(MCMC["alphas"]));
  field<mat> betas = List2Field_mat(as<List>(MCMC["betas"]));
  for (uword j = 0; j < betas.n_elem; ++j) betas.at(j) = trans(betas.at(j));
  mat sigmas = trans(as<mat>(MCMC["sigmas"]));
  cube D = as<cube>(MCMC["D"]);
  uword K = D.n_slices;
  uword q = D.slice(0).n_rows;
  mat sds(q, K);
  cube L(q, q, K);
  for (uword i = 0; i < K; ++i) {
    mat D_i = D.slice(i);
    sds.col(i) = sqrt(D_i.diag());
    mat R = cov2cor(D_i);
    L.slice(i) = chol(R);
  }
  //////////////////////////////
  // Sampling Random Effects //
  /////////////////////////////
  uword n_samples = as<uword>(control["n_samples"]);
  uword n_iter = as<uword>(control["n_iter"]);
  uword n_b = b_mat.n_rows;
  uword nRE = b_mat.n_cols;
  mat scale_b = mat(n_b,  b_mat.n_cols, fill::ones) * 0.2;
  //
  field<vec> betas_it(betas.n_elem);
  cube out(n_b, nRE, n_samples, fill::zeros);
  for (uword it = 0; it < n_samples; ++it) {
    vec bs_gammas_it = bs_gammas.col(it);
    vec gammas_it = gammas.col(it);
    vec alphas_it = alphas.col(it);
    for (uword i = 0; i < betas.n_elem; ++i) betas_it.at(i) = betas.at(i).col(it);
    vec sigmas_it = sigmas.col(it);
    mat L_it = L.slice(it);
    vec sds_it = sds.col(it);
    ///////////////////////
    vec W0H_bs_gammas = W0_H * bs_gammas_it;
    vec W0h_bs_gammas(W0_h.n_rows, fill::zeros);
    vec W0H2_bs_gammas(W0_H2.n_rows, fill::zeros);
    if (any_event) {
      W0h_bs_gammas = W0_h * bs_gammas_it;
    }
    if (any_interval) {
      W0H2_bs_gammas = W0_H2 * bs_gammas_it;
    }
    vec WH_gammas(W0_H.n_rows, fill::zeros);
    vec Wh_gammas(W0_h.n_rows, fill::zeros);
    vec WH2_gammas(W0_H2.n_rows);
    if (any_gammas) {
      WH_gammas = W_H * gammas_it;
    }
    if (any_gammas && any_event) {
      Wh_gammas = W_h * gammas_it;
    }
    if (any_gammas && any_interval) {
      WH2_gammas = W_H2 * gammas_it;
    }
    mat Wlong_H =
      calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, Wlong_sds, betas_it, b,
                      id_H_, FunForms, FunForms_ind, Funs_FunForms);
    vec WlongH_alphas = Wlong_H * alphas_it;
    mat Wlong_h(W0_h.n_rows, Wlong_H.n_cols, fill::zeros);
    vec Wlongh_alphas(W0_h.n_rows, fill::zeros);
    if (any_event) {
      Wlong_h =
        calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, Wlong_sds, betas_it, b,
                        id_h, FunForms, FunForms_ind, Funs_FunForms);
      Wlongh_alphas = Wlong_h * alphas_it;
    }
    mat Wlong_H2(W0_H2.n_rows, Wlong_H.n_cols, fill::zeros);
    vec WlongH2_alphas(W0_H2.n_rows, fill::zeros);
    if (any_interval) {
      Wlong_H2 =
        calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, Wlong_sds, betas_it,
                        b, id_H_, FunForms, FunForms_ind, Funs_FunForms);
      WlongH2_alphas = Wlong_H2 * alphas_it;
    }
    vec logLik_surv =
      log_surv_old(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas, //!! new //?? update later
               WH_gammas, Wh_gammas, WH2_gammas,
               WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
               log_Pwk, log_Pwk2, indFast_H, indFast_h,
               which_event, which_right_event, which_left,
               any_interval, which_interval);
    ///
    field<vec> eta = linpred_mixed(X, betas_it, Z, b, idL);
    vec logLik_long = log_long(y, eta, sigmas_it, extra_parms, families,
                               links, ids, unq_idL, n_b);
    ///
    vec logLik_re = log_re(b_mat, L_it, sds_it);
    // calculate the denominator
    vec denominator_b =
      logLik_long + logLik_surv + logLik_re;
    for (uword i = 0; i < n_iter; ++i) {
      for (uword j = 0; j < nRE; ++j) {
        //mat proposed_b_mat = propose_rnorm_mat(b_mat, scale_b, j);
        mat proposed_b_mat = b_mat;
        proposed_b_mat.col(j) = scale_b.col(j) % randn(b_mat.n_rows, 1) + b_mat.col(j);
        field<mat> proposed_b = mat2field(proposed_b_mat, ind_RE);
        //
        field<vec> eta_proposed =
          linpred_mixed(X, betas_it, Z, proposed_b, idL);
        vec logLik_long_proposed =
          log_long(y, eta_proposed, sigmas_it, extra_parms,
                   families, links, ids, unq_idL, n_b);
        //
        mat Wlong_H_proposed =
          calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, Wlong_sds,
                          betas_it, proposed_b, id_H_, FunForms,
                          FunForms_ind, Funs_FunForms);
        vec WlongH_alphas_proposed = Wlong_H_proposed * alphas_it;

        mat Wlong_h_proposed(Wlong_h.n_rows, Wlong_h.n_cols);
        vec Wlongh_alphas_proposed(Wlongh_alphas.n_rows);
        if (any_event) {
          Wlong_h_proposed =
            calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, Wlong_sds,
                            betas_it, proposed_b, id_h, FunForms,
                            FunForms_ind, Funs_FunForms);
          Wlongh_alphas_proposed = Wlong_h_proposed * alphas_it;
        }
        mat Wlong_H2_proposed(Wlong_H2.n_rows, Wlong_H2.n_cols);
        vec WlongH2_alphas_proposed(WlongH2_alphas.n_rows);
        if (any_interval) {
          Wlong_H2_proposed =
            calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, Wlong_sds, betas_it,
                            proposed_b, id_H_, FunForms, FunForms_ind,
                            Funs_FunForms);
          WlongH2_alphas_proposed = Wlong_H2_proposed * alphas_it;
        }
        //
        vec logLik_surv_proposed =
          log_surv_old(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas, //!! new //?? update later
                       WH_gammas, Wh_gammas, WH2_gammas,
                       WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                       log_Pwk, log_Pwk2, indFast_H, indFast_h,
                       which_event, which_right_event, which_left,
                       any_interval, which_interval);
        //
        vec logLik_re_proposed = log_re(proposed_b_mat, L_it, sds_it);
        //
        vec numerator_b =
          logLik_long_proposed + logLik_surv_proposed + logLik_re_proposed;
        vec log_ratio = numerator_b - denominator_b;
        for (uword k = 0; k < n_b; ++k) {
          double acc_k(0.0);
          if (std::isfinite(log_ratio.at(k)) &&
              exp(log_ratio.at(k)) > R::runif(0.0, 1.0)) {
            acc_k = 1.0;
            b_mat.row(k) = proposed_b_mat.row(k);
            denominator_b.at(k) = numerator_b.at(k);
            logLik_long.at(k) = logLik_long_proposed.at(k);
            logLik_surv.at(k) = logLik_surv_proposed.at(k);
            logLik_re.at(k) = logLik_re_proposed.at(k);
            uword first_H = GK_k * ni_event.at(k, 0);
            uword last_H = GK_k * ni_event.at(k, 1) - 1;
            Wlong_H.rows(first_H, last_H) = Wlong_H_proposed.rows(first_H, last_H);
            WlongH_alphas.rows(first_H, last_H) =
              WlongH_alphas_proposed.rows(first_H, last_H);
            if (any_event) {
              uword fitst_h = ni_event.at(k, 0);
              uword last_h = ni_event.at(k, 1) - 1;
              Wlong_h.rows(fitst_h, last_h) = Wlong_h_proposed.rows(fitst_h, last_h);
              Wlongh_alphas.rows(fitst_h, last_h) =
                Wlongh_alphas_proposed.rows(fitst_h, last_h);
            }
            if (any_interval) {
              Wlong_H2.rows(first_H, last_H) = Wlong_H2_proposed.rows(first_H, last_H);
              WlongH2_alphas.rows(first_H, last_H) =
                WlongH2_alphas_proposed.rows(first_H, last_H);
            }
          }
          if (i > 9) scale_b.at(k, j) = robbins_monro(scale_b.at(k, j), acc_k, i);
        }
      }
      b = mat2field(b_mat, ind_RE);
      eta = linpred_mixed(X, betas_it, Z, b, idL);
    }
    out.slice(it) = b_mat;
  }
  return out;
}

// [[Rcpp::export]]
arma::mat cum_haz (const List &Data, const List &MCMC) {
  ////////////////////////////
  // MCMC Sample Parameters //
  ///////////////////////////
  cube b = as<cube>(MCMC["b"]);
  mat bs_gammas = trans(as<mat>(MCMC["bs_gammas"]));
  mat gammas = trans(as<mat>(MCMC["gammas"]));
  mat alphas = trans(as<mat>(MCMC["alphas"]));
  field<mat> betas = List2Field_mat(as<List>(MCMC["betas"]));
  for (uword j = 0; j < betas.n_elem; ++j) betas.at(j) = trans(betas.at(j));
  vec W_std_gammas = as<vec>(MCMC["W_std_gammas"]);
  vec Wlong_std_alphas = as<vec>(MCMC["Wlong_std_alphas"]);
  uword n_samples = bs_gammas.n_cols;

  ////////////////////
  // index vectors //
  ///////////////////
  vec log_Pwk = as<vec>(Data["log_Pwk"]);
  uvec id_H = as<uvec>(Data["id_H"]) - 1;
  uvec id_h = as<uvec>(Data["id_h"]) - 1;
  uvec indFast_H = create_fast_ind(id_H + 1);
  uvec id_H_ = as<uvec>(Data["id_H_"]) - 1;
  uvec indFast_h = create_fast_ind(id_h + 1);
  field<uvec> ind_RE = List2Field_uvec(as<List>(Data["ind_RE"]), true);

  //////////////////////
  // Design matrices //
  /////////////////////
  mat W0_H = as<mat>(Data["W0_H"]);
  mat W_H = as<mat>(Data["W_H"]);
  field<mat> X_H = List2Field_mat(as<List>(Data["X_H"]));
  field<mat> Z_H = List2Field_mat(as<List>(Data["Z_H"]));
  field<mat> U_H = List2Field_mat(as<List>(Data["U_H"]));
  mat Wlong_bar = docall_cbindL(as<List>(Data["Wlong_bar"]));
  mat Wlong_sds = docall_cbindL(as<List>(Data["Wlong_sds"]));
  Wlong_bar = Wlong_bar.zeros();
  Wlong_sds = Wlong_sds.ones();

  bool any_gammas = as<bool>(Data["any_gammas"]);
  field<uvec> FunForms = List2Field_uvec(as<List>(Data["FunForms_cpp"]), true);
  field<uvec> FunForms_ind = List2Field_uvec(as<List>(Data["FunForms_ind"]), true);
  List Funs_FunForms = as<List>(Data["Funs_FunForms"]);

  ///////////////////////
  // Calculation CIFs //
  //////////////////////
  field<vec> betas_it(betas.n_elem);
  uvec H_unq = unique(id_H);
  mat out(H_unq.n_rows, n_samples, fill::zeros);
  for (uword it = 0; it < n_samples; ++it) {
    field<mat> b_it = mat2field(b.slice(it), ind_RE);
    vec bs_gammas_it = bs_gammas.col(it);
    vec gammas_it = gammas.col(it);
    vec alphas_it = alphas.col(it);
    for (uword i = 0; i < betas.n_elem; ++i) betas_it.at(i) = betas.at(i).col(it);
    ///////////////////////
    vec W0H_bs_gammas = (W0_H * bs_gammas_it) - W_std_gammas.at(it) - Wlong_std_alphas.at(it);
    vec WH_gammas(W0_H.n_rows, fill::zeros);
    if (any_gammas) {
      WH_gammas = W_H * gammas_it;
    }
    mat Wlong_H =
      calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, Wlong_sds, betas_it, b_it,
                      id_H_, FunForms, FunForms_ind, Funs_FunForms);
    vec WlongH_alphas = Wlong_H * alphas_it;
    ///////////////////////
    vec lambda_H = W0H_bs_gammas + WH_gammas + WlongH_alphas;
    out.col(it) = group_sum(exp(log_Pwk + lambda_H), indFast_H);
  }
  return out;
}



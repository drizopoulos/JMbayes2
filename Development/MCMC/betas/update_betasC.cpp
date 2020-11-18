#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

vec docall_rbindF (const field<vec> &F) { // binds a field of vectors into one vector
  uword n = F.n_elem;
  uword nrows = 0;
  uvec rows(n);
  for (uword i = 0; i < n; i++) {
    rows.at(i) = F.at(i).n_rows;
    nrows += rows.at(i);
  }
  vec V(nrows);
  uword ii = 0;
  for (uword i = 0; i < n; i++) {
    V.rows(ii, ii - 1 + rows.at(i)) = F.at(i);
    ii += rows.at(i);
  }
  return V;
}

// ?? The function below could be optimized to update in blocks, rather then element by element
mat add_zero_colrows (const mat &M, // adds zero-rows and/or zero-cols to a matrix M
                      const uword &nrows, // n_rows in the target matrix
                      const uword &ncols, // n_cols in the target matrix
                      const uvec &rows_ind,    // ind where to place the M's rows (zero-rows will be 'added' to the absent ind). the number of ind must match the M's n_rows
                      const uvec &cols_ind) { // ind where to place the M's cols (zero-cols will be 'added' to the absent ind). the number of ind must match the M's n_cols
  
  mat Res(nrows, ncols, fill::zeros);
  uword M_nrows = M.n_rows;
  uword M_ncols = M.n_cols;
  
  for(uword i = 0; i < M_nrows; i++) { // by row
    for(uword j = 0; j < M_ncols; j++) { // by col
      Res.at(rows_ind.at(i), cols_ind.at(j)) = M.at(i, j);
    }
  }
  return Res;
}

mat add_zero_rows (const mat &M, // adds zero-rows to a matrix M
                   const uword &nrows, // n_rows in the target matrix
                   const uvec &rows_ind) { // ind where to place the M's rows (zero-rows will be 'added' to the absent ind). the number of ind must match the M's n_rows
  
  mat Res(nrows, M.n_cols, fill::zeros);
  uword M_nrows = M.n_rows;
  
  for(uword j = 0; j < M_nrows; j++) {
    Res.row(rows_ind.at(j)) = M.row(j);
  }
  return Res;
}

mat rank1_update (const mat &U, // performs rank-1 update. If U = chol(M), returns chol(M + v * v.t())
                  const vec &v) {
  
  uword n = v.n_elem;
  mat Res = U;
  vec v2  = v;
  
  for (uword i = 0; i < n; i++) {
    double r = pow( pow(Res.at(i, i), 2) + pow(v2.at(i), 2), 0.5);
    double c = r / Res.at(i, i);
    double s = v2.at(i) / Res.at(i, i);
    Res.at(i, i) = r;
    
    if (i < n-1) {
      Res.submat(i, i + 1, i, n - 1) = (Res.submat(i, i + 1, i, n - 1) + s * v2.rows(i + 1, n - 1).t()) / c;
      v2.rows(i + 1, n - 1) = c * v2.rows(i + 1, n - 1) - s * Res.submat(i, i + 1, i, n - 1).t();
    }
  }
  return Res;
}

mat chol_update(const mat &U, // If U = chol(M), returns chol(M.submat(keep, keep))
                const uvec &keep) { // keep must be a sorted vector, i.e, {2, 4, 5}, and counts from 0
  
  // to improve: later we can try to extend this approach further to obtain inv(U) from the required inv(U_i)
  
  uvec rem = regspace<uvec>(0,  U.n_cols - 1); rem.shed_rows(keep); // cols-rows to remove
  mat Res = U;
  uword n = rem.n_elem;
  
  for (uword i = 0; i < n; i++) { // rank-1 update for each col-row to be removed
    
    uword last_col = Res.n_cols - 1;
    
    if(rem.at(i) < last_col) {
      Res.submat(rem.at(i) + 1, rem.at(i) + 1, last_col, last_col) = rank1_update(Res.submat(rem.at(i) + 1, rem.at(i) + 1, last_col, last_col),
                 Res.submat(rem.at(i), rem.at(i) + 1, rem.at(i), last_col).t());
    }
    
    Res.shed_row(rem.at(i));
    Res.shed_col(rem.at(i));
    rem = rem - 1;
  }
  return Res;
}

mat rnorm_mat (const uword& rows, const uword& cols) {
  mat out(rows, cols);
  out.each_col([&](vec& x) {x = as<vec>(rnorm(rows)); } );
  return out;
}

mat propose_mvnorm_vec (const int &n, // number of samples
                        const mat &U, // upper triangular Choleski factorization of the vcov matrix
                        const double &scale) {
  uword ncols = U.n_cols;
  mat res(n, ncols);
  res = scale * (rnorm_mat(n, ncols) * U);
  return res.t(); // j-th column reports the j-th sample
}

field<vec> vec2field (const vec &betas, const field<uvec> &ind_FE) {
  uword n = ind_FE.n_elem;
  field<vec> out(n);
  for (uword i = 0; i < n; i++) {
    out.at(i) = betas.rows(ind_FE.at(i));
  }
  return out;
}

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

mat docall_cbindF (const field<mat> &Mats) {
  uword n = Mats.n_elem;
  uvec ncols(n);
  for (uword k = 0; k < n; ++k) {
    ncols.at(k) = Mats.at(k).n_cols;
  }
  uword N = sum(ncols);
  uword col_start = 0;
  uword col_end = ncols.at(0) - 1;
  mat out(Mats.at(0).n_rows, N);
  for (uword k = 0; k < n; ++k) {
    if (k > 0) {
      col_start += ncols.at(k - 1);
      col_end += ncols.at(k);
    }
    out.cols(col_start, col_end) = Mats.at(k);
  }
  return out;
}

mat cov2cor (const mat &V) {
  vec Is = sqrt(1.0 / V.diag());
  mat out = V.each_col() % Is;
  out.each_row() %= Is.t();
  return out;
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

////////////////////////////////////////////////////////////////////////////////

void update_betas (field<vec> &betas, mat &res_betas, //mat &acceptance_betas,
                   //vec &scale_betas, field<vec> &eta, vec &logLik_long,
                   //vec &logLik_surv, mat &Wlong_H, mat &Wlong_h, mat &Wlong_H2,
                   //vec &WlongH_alphas, vec &Wlongh_alphas, vec &WlongH2_alphas,
                   const vec &prior_mean_betas_HC, const mat &prior_Tau_betas_HC,
                   const mat &b_mat, const mat &L, const vec &sds, const mat &X_dot,
                   const field<uvec> &ind_FE, // indices for the FE in res_betas[it,] belonging to the field betas. E.g., {{1,2,3}, {4, 5}, {6}}
                   const uvec &ind_FE_HC, // indices for the FE present in the HC (cols in res_betas)
                   const uvec &id_patt, // vector with the ids' outcome missing pattern
                   const field<uvec> &ind_RE_patt, // indices for the RE present in each outcome missing pattern (cols in D)
                   const field<uvec> &ind_FE_patt, // indices for the FE (in HC) present in each outcome missi
                   const uword &it,/*
                   const uvec &has_tilde_betas,
                   const field<mat> &X,
                   const field<mat> &Z,
                   const field<mat> &b,
                   const field<uvec> &idL,*/
                   const field<mat> &y/*
                   const vec &sigmas,
                   const vec &extra_parms,
                   const CharacterVector &families,
                   const CharacterVector &links,
                   const field<uvec> &idL_lp_fast,
                   const field<vec> &prior_mean_betas_nHC,
                   field<mat> &prior_Tau_betas_nHC,
                   const field<mat> &chol_vcov_prop_betas,
                   const field<uvec> &x_notin_z,
                   const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
                   const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
                   const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
                   const mat &Wlong_bar,
                   const uvec &id_H_, const uvec &id_h,
                   const field<uvec> &FunForms,
                   const field<uvec> &FunForms_ind,
                   const vec &alphas,
                   const bool &any_event, const bool &any_interval,
                   const vec &W0H_bs_gammas, const vec &W0h_bs_gammas, const vec &W0H2_bs_gammas,
                   const vec &WH_gammas, const vec &Wh_gammas, const vec &WH2_gammas,
                   const vec &log_Pwk, const vec &log_Pwk2,
                   const uvec &id_H_fast, const uvec &id_h_fast,
                   const uvec &which_event, const uvec &which_right_event, const uvec &which_left,
                   const uvec &which_interval,
                   const field<uvec> &unq_idL*/
                   ) {
  /*
   * To improve: initialize variables outside loops, and check if it gets faster (applies for both Gibbs and MH)
   */
  
  uword n_b = b_mat.n_rows;
  uword n_outcomes = y.n_elem;
  
  // FE in HC - Gibbs sampling
  vec betas_vec = docall_rbindF(betas);
  uword patt_count = id_patt.max() + 1; // number of unique outcome-missing patterns
  uword p_HC = ind_FE_HC.n_elem; // number of HC-FE
  uword q = b_mat.n_cols; // number of RE
  mat sum_JXDXJ(p_HC, p_HC, fill::zeros); // sum for the posterior's parameters
  vec sum_JXDu(p_HC, fill::zeros); // sum for the posterior's parameters
  mat U = L.each_row() % sds.t(); // RE vcov matrix Cholesky factorization (upper)
  field<mat> D_inv(patt_count); // all unique vcov_inv matrices across the missing outcome patterns
  for (uword i = 0; i < n_b; ++i) { // i-th patient
    /* I'm assuming that n_b in the total number of patients, including those
     * who only have survival outcome. I.e., b_mat has rows for them too.
     */
    if (i < patt_count && !ind_RE_patt.at(i).is_empty()) {
      /* obtain all unique vcov_inv matrices required for the sums in the posterior parameters
       * & jumps the pattern in which the patient misses all longitudinal outcomes
       */
      mat U_patt_inv = inv(trimatu(chol_update(U, ind_RE_patt.at(i))));
      D_inv.at(i) =  U_patt_inv * U_patt_inv.t(); // mat
    }
    
    
    uword patt_i = id_patt.at(i); // id missing outcome pattern
    if (ind_FE_patt.at(patt_i).is_empty()) {continue;} // skip ids without longitudinal outcomes
    uvec ind_FE_i = ind_FE_patt.at(patt_i);
    uvec ind_RE_i = ind_RE_patt.at(patt_i);
    mat X_dot_i = X_dot.rows(i * q, (i + 1) * q - 1);
    X_dot_i = X_dot_i.submat(ind_RE_i, ind_FE_i);
    vec b_i = b_mat.row(i).t();
    vec u_i = b_i.rows(ind_RE_i) + X_dot_i * betas_vec.rows(ind_FE_i);
    mat D_inv_i = D_inv.at(patt_i);
    mat XD_i = X_dot_i.t() * D_inv_i;
    mat XDX_i = XD_i * X_dot_i;
    sum_JXDu  += add_zero_rows(XD_i*u_i, p_HC, ind_FE_i);
    sum_JXDXJ += add_zero_colrows(XDX_i, p_HC, p_HC, ind_FE_i, ind_FE_i);
  }
  mat Sigma_1 = inv(prior_Tau_betas_HC + sum_JXDXJ); // improve via Cholesky decomposition 
  vec mean_1 = Sigma_1 * (prior_Tau_betas_HC * prior_mean_betas_HC + sum_JXDu);
  mat U_1 = chol(Sigma_1);
  betas_vec.rows(ind_FE_HC) = propose_mvnorm_vec(1, U_1, 1.0) + mean_1;
  betas = vec2field(betas_vec, ind_FE);
  
  
  Rcout << "mean_1: \n"  << mean_1 << "\n";
  Rcout << "Sigma_1: \n\n" << Sigma_1 << "\n";
  Rcout << "D_inv.at(0): " << D_inv.at(0) << "\n";
  Rcout << "D_inv.at(1): \n" << D_inv.at(1) << "\n";
  Rcout << "D_inv.at(2): \n" << D_inv.at(2) << "\n";
  Rcout << "D_inv.at(3): \n" << D_inv.at(3) << "\n";
  
  /*
  // update eta
  eta = linpred_mixed(X, betas, Z, b, idL);
  // update logLik_surv
  Wlong_H =
    calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, betas, b, id_H_, FunForms,
                    FunForms_ind);
  WlongH_alphas = Wlong_H * alphas;
  if (any_event) {
    Wlong_h =
      calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, betas, b, id_h, FunForms,
                      FunForms_ind);
    Wlongh_alphas = Wlong_h * alphas;
  }
  if (any_interval) {
    Wlong_H2 =
      calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, betas, b, id_H_, FunForms,
                      FunForms_ind);
    WlongH2_alphas = Wlong_H2 * alphas;
  }
  logLik_surv = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                         WH_gammas, Wh_gammas, WH2_gammas,
                         WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                         log_Pwk, log_Pwk2, id_H_fast, id_h_fast,
                         which_event, which_right_event, which_left,
                         any_interval, which_interval);
  
  ///////////////////////////////////////////////////////////////////////////////
  // FE outside HC - Metropolis-Hastings sampling
  if (any(has_tilde_betas)) {
    for (uword j = 0; j < n_outcomes; ++j) { // j-th outcome
      if (!has_tilde_betas.at(j)) {continue;} // skip outcomes without nHC-FE
      uvec ind_j = x_notin_z.at(j);
      // denominator
      double sum_logLik_long_j =
        sum(log_long_i(y.at(j), eta.at(j), sigmas.at(j), extra_parms.at(j),
                       std::string(families[j]), std::string(links[j]),
                       idL_lp_fast.at(j)));
      vec ll(prior_mean_betas_nHC.at(j).n_rows, fill::ones);
      double prior =
        logPrior(betas.at(j).rows(ind_j), prior_mean_betas_nHC.at(j),
                 prior_Tau_betas_nHC.at(j), ll, 1.0, false);
      double denominator_j = sum_logLik_long_j + sum(logLik_surv) + prior;
      // proposal
      field<vec> betas_prop = betas;
      betas_prop.at(j).rows(ind_j) += propose_mvnorm_vec(1, chol_vcov_prop_betas.at(j),
                    scale_betas.at(j));
      double prior_prop =
        logPrior(betas_prop.at(j).rows(ind_j), prior_mean_betas_nHC.at(j),
                 prior_Tau_betas_nHC.at(j), ll, 1.0, false);
      // logLik_long proposal
      field<vec> eta_prop = linpred_mixed(X, betas_prop, Z, b, idL);
      double sum_logLik_long_j_prop =
        sum(log_long_i(y.at(j), eta_prop.at(j), sigmas.at(j), extra_parms.at(j),
                       std::string(families[j]), std::string(links[j]),
                       idL_lp_fast.at(j)));
      // logLik_surv proposal
      mat Wlong_H_prop =
        calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, betas_prop, b, id_H_,
                        FunForms, FunForms_ind);
      vec WlongH_alphas_prop = Wlong_H_prop * alphas;
      mat Wlong_h_prop;
      vec Wlongh_alphas_prop;
      if (any_event) {
        Wlong_h_prop =
          calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, betas_prop, b, id_h,
                          FunForms, FunForms_ind);
        Wlongh_alphas_prop = Wlong_h_prop * alphas;
      }
      mat Wlong_H2_prop;
      vec WlongH2_alphas_prop;
      if (any_interval) {
        Wlong_H2_prop =
          calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, betas_prop, b, id_H_,
                          FunForms, FunForms_ind);
        WlongH2_alphas_prop = Wlong_H2_prop * alphas;
      }
      vec logLik_surv_prop =
        log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                 WH_gammas, Wh_gammas, WH2_gammas,
                 WlongH_alphas_prop, Wlongh_alphas_prop, WlongH2_alphas_prop,
                 log_Pwk, log_Pwk2, id_H_fast, id_h_fast,
                 which_event, which_right_event, which_left,
                 any_interval, which_interval);
      // numerator
      double numerator_j =
        sum_logLik_long_j_prop + sum(logLik_surv_prop) + prior_prop;
      // Hastings ratio
      double log_ratio_j = numerator_j - denominator_j;
      if(std::isfinite(log_ratio_j) &&
         std::exp(log_ratio_j) > R::runif(0.0, 1.0)) {
        acceptance_betas.at(it, j) = 1;
        betas.at(j) = betas_prop.at(j);
        eta = eta_prop;
        Wlong_H = Wlong_H_prop;
        Wlong_h = Wlong_h_prop;
        Wlong_H2 = Wlong_H2_prop;
        WlongH_alphas = WlongH_alphas_prop;
        Wlongh_alphas = Wlongh_alphas_prop;
        WlongH2_alphas = WlongH2_alphas_prop;
        logLik_surv = logLik_surv_prop;
        if(it > 19) {
          scale_betas.at(j) =
            robbins_monro(scale_betas.at(j), acceptance_betas.at(it, j), it,
                          0.25);
        }
      }
    }
  }
  // update logLik_long with all betas
  logLik_long = log_long(y, eta, sigmas, extra_parms,
                         families, links, idL_lp_fast, unq_idL, n_b);
  
  */
  // save all results
  res_betas.row(it) = docall_rbindF(betas).t();
}

// [[Rcpp::export]]
List update_betas_Gibbs (List control, List initial_values, List priors, List model_data) {
  
  // control
  uword n_iter = as<uword>(control["n_iter"]);
  
  // initial values
  field<vec> betas = List2Field_vec(as<List>(initial_values["betas"]));
  field<mat> b = List2Field_mat(as<List>(initial_values["b"]));
  mat D = as<mat>(initial_values["D"]);
  
  //
  vec betas_vec = docall_rbindF(betas);
  uword n_betas = betas_vec.n_rows;
  mat res_betas(n_iter, n_betas, fill::zeros);
  mat b_mat = docall_cbindF(b);
  uword n_b = b_mat.n_rows;
  vec sds = sqrt(D.diag());
  mat R = cov2cor(D);
  mat L = chol(R);
  
  // priors
  vec prior_mean_betas_HC = as<vec>(priors["mean_betas_HC"]);
  mat prior_Tau_betas_HC = as<mat>(priors["Tau_betas_HC"]);
  
  //model_data
  mat X_dot = as<mat>(model_data["X_dot"]);
  field<uvec> ind_FE = List2Field_uvec(as<List>(model_data["ind_FE"]), true);
  uvec ind_FE_HC = as<uvec>(model_data["ind_FE_HC"]) - 1;
  uvec id_patt = as<uvec>(model_data["id_patt"]) - 1;
  field<uvec> ind_RE_patt = List2Field_uvec(as<List>(model_data["ind_RE_patt"]), true);
  field<uvec> ind_FE_patt = List2Field_uvec(as<List>(model_data["ind_FE_patt"]), true);
  field<mat> y = List2Field_mat(as<List>(model_data["y"]));

  //////////////////////////////////////////////////////////////////////////////

  for (uword it = 0; it < n_iter; ++it) {
    
    update_betas (betas,res_betas, 
                  prior_mean_betas_HC, prior_Tau_betas_HC,
                  b_mat, L, sds, X_dot,
                  ind_FE, // indices for the FE in res_betas[it,] belonging to the field betas. E.g., {{1,2,3}, {4, 5}, {6}}
                  ind_FE_HC, // indices for the FE present in the HC (cols in res_betas)
                  id_patt, // vector with the ids' outcome missing pattern
                  ind_RE_patt, // indices for the RE present in each outcome missing pattern (cols in D)
                  ind_FE_patt, // indices for the FE (in HC) present in each outcome missi
                  it,
                  y);
  }
  
  return List::create(
  );
}


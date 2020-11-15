#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

field<vec> linpred_mixed (const field<mat> &X, const field<vec> &betas,
                          const field<mat> &Z, const field<mat> &b,
                          const field<uvec> &id) {
  uword n_outcomes = X.n_elem;
  field<vec> out(n_outcomes);
  for (uword i = 0; i < n_outcomes; ++i) {
    mat X_i = X.at(i);
    vec betas_i = betas.at(i);
    mat Z_i = Z.at(i);
    mat b_i = b.at(i);
    uvec id_i = id.at(i);
    out.at(i) = X_i * betas_i + arma::sum(Z_i % b_i.rows(id_i), 1);
  }
  return out;
}


vec mu_fun (const vec &eta, const std::string &link) {
  uword n = eta.n_rows;
  vec out(n);
  if (link == "identity") {
    out = eta;
  } else if (link == "inverse") {
    out = 1.0 / eta;
  } else if (link == "logit") {
    out = 1.0 / (1.0 + trunc_exp(- eta));
  } else if (link == "probit") {
    out = normcdf(eta);
  } else if (link == "cloglog") {
    out = - trunc_exp(- trunc_exp(eta)) + 1.0;
  } else if (link == "log") {
    out = trunc_exp(eta);
  }
  return out;
}

vec log_dnorm (const vec &x, const vec &mu, const double &sigma) {
  vec sigmas(x.n_rows);
  sigmas.fill(sigma);
  vec out = log_normpdf(x, mu, sigmas);
  return out;
}

vec log_pnorm (const vec &x, const vec &mu, const double &sigma,
               const int lower_tail = 1) {
  uword n = x.n_rows;
  vec out(n);
  for (uword i = 0; i < n; ++i) {
    out.at(i) = R::pnorm(x.at(i), mu.at(i), sigma, lower_tail, 1);
  }
  return out;
}

vec log_dt (const vec &x, const double &df) {
  uword n = x.n_rows;
  vec out(n);
  for (uword i = 0; i < n; ++i) {
    out.at(i) = R::dt(x.at(i), df, 1);
  }
  return out;
}

vec log_dgamma (const vec &x, const double &shape, const vec &scale) {
  uword n = x.n_rows;
  vec out(n);
  for (uword i = 0; i < n; ++i) {
    out.at(i) = R::dgamma(x.at(i), shape, scale.at(i), 1);
  }
  return out;
}

vec log_dbeta (const vec &x, const vec &shape1, const vec &shape2) {
  uword n = x.n_rows;
  vec out(n);
  for (uword i = 0; i < n; ++i) {
    out.at(i) = R::dbeta(x.at(i), shape1.at(i), shape2.at(i), 1);
  }
  return out;
}

vec log_dbinom (const vec &x, const vec &size, const vec &prob) {
  uword n = x.n_rows;
  vec out(n);
  for (uword i = 0; i < n; ++i) {
    out.at(i) = R::dbinom(x.at(i), size.at(i), prob.at(i), 1);
  }
  return out;
}

vec log_dpois (const vec &x, const vec &lambda) {
  uword n = x.n_rows;
  vec out(n);
  for (uword i = 0; i < n; ++i) {
    out.at(i) = R::dpois(x.at(i), lambda.at(i), 1);
  }
  return out;
}


vec lbeta_arma (const vec &a, const vec &b) {
  uword n = a.n_rows;
  vec out(n);
  for (uword i = 0; i < n; ++i) {
    out.at(i) = R::lbeta(a.at(i), b.at(i));
  }
  return out;
}

vec lchoose_arma (const vec &n, const vec &k) {
  uword n_ = n.n_rows;
  vec out(n_);
  for (uword i = 0; i < n_; ++i) {
    out.at(i) = R::lchoose(n.at(i), k.at(i));
  }
  return out;
}

vec log_dbbinom (const vec &x, const vec &size, const vec &prob,
                 const double &phi) {
  vec A = phi * prob;
  vec B = phi * (1.0 - prob);
  vec log_numerator = lbeta_arma(x + A, size - x + B);
  vec log_denominator = lbeta_arma(A, B);
  vec fact = lchoose_arma(size, x);
  vec out = fact + log_numerator - log_denominator;
  return out;
}

vec log_dnbinom (const vec &x, const vec &mu, const double &size) {
  vec log_mu_size = log(mu + size);
  vec comp1 = lgamma(x + size) - lgamma(size) - lgamma(x + 1);
  vec comp2 = size * log(size) - size * log_mu_size;
  vec comp3 = x % (log(mu) - log_mu_size);
  vec out = comp1 + comp2 + comp3;
  return out;
}

vec group_sum (const vec &x, const uvec &ind) {
  vec cumsum_x = cumsum(x);
  vec out = cumsum_x.rows(ind);
  out.insert_rows(0, 1);
  out = diff(out);
  return out;
}

vec log_long_i (const mat &y_i, const vec &eta_i, const double &sigma_i,
                const double &extr_prm_i, const std::string &fam_i,
                const std::string &link_i, const uvec &idFast_i) {
  uword N = y_i.n_rows;
  vec log_contr(N);
  vec mu_i = mu_fun(eta_i, link_i);
  if (fam_i == "gaussian") {
    log_contr = log_dnorm(y_i, mu_i, sigma_i);
  } else if (fam_i == "Student-t") {
    log_contr = log_dt((y_i - mu_i) / sigma_i, extr_prm_i) - std::log(sigma_i);
  } else if (fam_i == "beta") {
    log_contr = log_dbeta(y_i, mu_i * sigma_i, sigma_i * (1.0 - mu_i));
  } else if (fam_i == "Gamma") {
    log_contr = log_dgamma(y_i, sigma_i, mu_i / sigma_i);
  } else if (fam_i == "unit Lindley") {
    vec theta = 1.0 / mu_i - 1.0;
    vec comp1 = 2.0 * log(theta) - log(1.0 + theta);
    vec comp2 = - 3.0 * log(1.0 - y_i);
    vec comp3 = - (theta * y_i) / (1.0 - y_i);
    log_contr = comp1 + comp2 + comp3;
  } else if (fam_i == "censored normal") {
    uvec ind0 = find(y_i.col(1) == 0);
    uvec ind1 = find(y_i.col(1) == 1);
    uvec ind2 = find(y_i.col(1) == 2);
    vec yy = y_i.col(0);
    log_contr.rows(ind0) = log_dnorm(yy.rows(ind0), mu_i.rows(ind0), sigma_i);
    log_contr.rows(ind1) = log_pnorm(yy.rows(ind1), mu_i.rows(ind1), sigma_i);
    log_contr.rows(ind2) = log_pnorm(yy.rows(ind2), mu_i.rows(ind2), sigma_i, 0);
  } else if (fam_i == "binomial") {
    uword k = y_i.n_cols;
    if (k == 2) {
      // y_i.col(1) has been set to y_i.col(0) + y_i.col(1)
      // in jm_fit(), i.e., y_i.col(1) is already the number of trials
      // not the number of failures
      log_contr = log_dbinom(y_i.col(0), y_i.col(1), mu_i);
    } else {
      log_contr = y_i % log(mu_i) + (1.0 - y_i) % log(1.0 - mu_i);
    }
  } else if (fam_i == "poisson") {
    log_contr = log_dpois(y_i, mu_i);
  } else if (fam_i == "negative binomial") {
    log_contr = log_dnbinom(y_i, mu_i, sigma_i);
  }  else if (fam_i == "beta binomial") {
    uword k = y_i.n_cols;
    if (k == 2) {
      // y_i.col(1) has been set to y_i.col(0) + y_i.col(1)
      // in jm_fit(), i.e., y_i.col(1) is already the number of trials
      // not the number of failures
      log_contr = log_dbbinom(y_i.col(0), y_i.col(1), mu_i, sigma_i);
    } else {
      vec ones(N, fill::ones);
      log_contr = log_dbbinom(y_i, ones, mu_i, sigma_i);
    }
  }
  vec out = group_sum(log_contr, idFast_i);
  return out;
}

static double const log2pi = std::log(2.0 * M_PI);

void inplace_UpperTrimat_mult (rowvec &x, const mat &trimat) {
  // in-place multiplication of x with an upper triangular matrix trimat
  // because in-place assignment is much faster but careful in use because
  // it changes the input vector x, i.e., not const
  uword const n = trimat.n_cols;
  for (uword j = n; j-- > 0;) {
    double tmp(0.0);
    for (uword i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x.at(i);
    x.at(j) = tmp;
  }
}

vec log_dmvnrm_chol (const mat &x, const mat &L) {
  // fast log density of the multivariate normal distribution
  // L is the Cholesky factor of the covariance matrix.
  using arma::uword;
  uword const n = x.n_rows, xdim = x.n_cols;
  vec out(n);
  mat V = inv(trimatu(L));
  double const log_det = sum(log(V.diag())),
    constants = -(double)xdim / 2.0 * log2pi,
    other_terms = constants + log_det;
  rowvec z_i(xdim);
  for (uword i = 0; i < n; i++) {
    z_i = x.row(i);
    inplace_UpperTrimat_mult(z_i, V);
    out.at(i) = other_terms - 0.5 * dot(z_i, z_i);
  }
  return out;
}

field<mat> linpred_surv (const field<mat> &X, const field<vec> &betas,
                         const field<mat> &Z, const field<mat> &b,
                         const uvec &id) {
  uword n_outcomes = X.n_elem;
  field<mat> out(n_outcomes);
  for (uword i = 0; i < n_outcomes; ++i) {
    mat X_i = X.at(i);
    vec betas_i = betas.at(i);
    mat Z_i = Z.at(i);
    mat b_i = b.at(i);
    uword n_betas = betas_i.n_rows;
    uword n_REs = b_i.n_cols;
    uword n_forms = X_i.n_cols / n_betas;
    mat out_i(X_i.n_rows, n_forms);
    out.at(i) = out_i;
    for (uword j = 0; j < n_forms; ++j) {
      mat X_ij = X_i.cols(j * n_betas, (j + 1) * n_betas - 1);
      mat Z_ij = Z_i.cols(j * n_REs, (j + 1) * n_REs - 1);
      out.at(i).col(j) = X_ij * betas_i +
        arma::sum(Z_ij % b_i.rows(id), 1);
    }
  }
  return out;
}

field<mat> create_Wlong (const field<mat> &eta, const field<uvec> &FunForms,
                         const field<mat> &U, const field<uvec> &ind) {
  uword n_outcomes = eta.n_elem;
  field<mat> out(n_outcomes);
  for (uword i = 0; i < n_outcomes; ++i) {
    mat eta_i = eta.at(i);
    uvec FF_i = FunForms.at(i);
    mat U_i = U.at(i);
    uvec ind_i = ind.at(i);
    mat Wlong_i(eta_i.n_rows, U_i.n_cols, fill::ones);
    Wlong_i.cols(FF_i) %= eta_i.cols(ind_i);
    out.at(i) = U_i % Wlong_i;
  }
  return out;
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

mat calculate_Wlong (const field<mat> &X, const field<mat> &Z,
                     const field<mat> &U, const mat &Wlong_bar,
                     const field<vec> &betas, const field<mat> &b,
                     const uvec &id, const field<uvec> &FunForms,
                     const field<uvec> &FunForms_ind) {
  field<mat> eta = linpred_surv(X, betas, Z, b, id);
  mat Wlong = docall_cbindF(create_Wlong(eta, FunForms, U, FunForms_ind));
  Wlong.each_row() -= Wlong_bar;
  return Wlong;
}

vec log_surv (const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
              const vec &W0H2_bs_gammas, const vec &WH_gammas,
              const vec &Wh_gammas, const vec &WH2_gammas,
              const vec &WlongH_alphas, const vec &Wlongh_alphas,
              const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
              const uvec &indFast_H, const uvec &indFast_h, const uvec &which_event,
              const uvec &which_right_event, const uvec &which_left,
              const bool &any_interval, const uvec &which_interval) {
  vec lambda_H = W0H_bs_gammas + WH_gammas + WlongH_alphas;
  vec H = group_sum(exp(log_Pwk + lambda_H), indFast_H);
  uword n = H.n_rows;
  vec lambda_h(n);
  lambda_h.rows(which_event) = W0h_bs_gammas.rows(which_event) +
    Wh_gammas.rows(which_event) + Wlongh_alphas.rows(which_event);
  vec out(n);
  out.rows(which_right_event) = - H.rows(which_right_event);
  out.rows(which_event) += lambda_h.rows(which_event);
  out.rows(which_left) = log1p(- exp(- H.rows(which_left)));
  vec lambda_H2(lambda_H.n_rows);
  vec H2(n);
  if (any_interval) {
    lambda_H2 = W0H2_bs_gammas + WH2_gammas + WlongH2_alphas;
    H2 = group_sum(exp(log_Pwk2 + lambda_H2), indFast_H);
    out.rows(which_interval) = - H.rows(which_interval) +
      log(- expm1(- H2.rows(which_interval)));
  }
  out = group_sum(out, indFast_h);
  return out;
}

vec log_re (const mat &b, const mat &L, const vec &sds) {
  mat chol_Sigma = L.each_row() % sds.t();
  vec out = log_dmvnrm_chol(b, chol_Sigma);
  return out;
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

double robbins_monro (const double &scale, const double &acceptance_it,
                      const int &it, const double &target_acceptance = 0.45) {
  double step_length = scale / (target_acceptance * (1.0 - target_acceptance));
  double out;
  if (acceptance_it > 0) {
    out = scale + step_length * (1 - target_acceptance) / (double)it;
  } else {
    out = scale - step_length * target_acceptance / (double)it;
  }
  return out;
}

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

vec log_long (const field<mat> &y, const field<vec> &eta, const vec &sigmas,
              const vec &extra_parms, const CharacterVector &families,
              const CharacterVector &links, const field<uvec> &idFast,
              const field<uvec> &unq_ids, const uword &n) {
  uword n_outcomes = y.size();
  vec out(n, fill::zeros);
  for (uword i = 0; i < n_outcomes; ++i) {
    mat y_i = y.at(i);
    vec eta_i = eta.at(i);
    double sigma_i = sigmas.at(i);
    double extr_prm_i = extra_parms.at(i);
    std::string fam_i = std::string(families[i]);
    std::string link_i = std::string(links[i]);
    uvec idFast_i = idFast.at(i);
    uvec unq_id_i = unq_ids.at(i);
    vec log_contr_i = log_long_i(y_i, eta_i, sigma_i, extr_prm_i, fam_i,
                                 link_i, idFast_i);
    out.rows(unq_id_i) += log_contr_i;
  }
  return out;
}

vec create_init_scale(const uword &n, const double &fill_val = 0.1) {
  vec out(n);
  out.fill(fill_val);
  return out;
}

field<vec> List2Field_vec (const List &Vecs) {
  uword n_list = Vecs.size();
  field<vec> res(n_list);
  for (uword i = 0; i < n_list; ++i) {
    res.at(i) = as<vec>(Vecs[i]);
  }
  return res;
}

field<mat> List2Field_mat (const List &Mats) {
  uword n_list = Mats.size();
  field<mat> res(n_list);
  for (uword i = 0; i < n_list; ++i) {
    res.at(i) = as<mat>(Mats[i]);
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

uvec create_fast_ind (const uvec &group) {
  uword l = group.n_rows;
  uvec ind = find(group.rows(1, l - 1) != group.rows(0, l - 2));
  uword k = ind.n_rows;
  ind.insert_rows(k, 1);
  ind.at(k) = l - 1;
  return ind;
}

mat docall_cbindL (const List &Mats_) {
  field<mat> Mats = List2Field_mat(Mats_);
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

////////////////////////////////////////////////////////////////////////////////

void update_betas (field<vec> &betas,
                   mat &res_betas,
                   mat &acceptance_betas,
                   vec &scale_betas,
                   field<vec> &eta,
                   vec &logLik_long,
                   vec &logLik_surv,
                   mat &Wlong_H, mat &Wlong_h, mat &Wlong_H2,
                   vec &WlongH_alphas, vec &Wlongh_alphas, vec &WlongH2_alphas,
                   const uword &it,
                   const uvec &has_tilde_betas,
                   const field<mat> &X, 
                   const field<mat> &Z,
                   const field<mat> &b,
                   const field<uvec> &idL,
                   const field<mat> &y,
                   const vec &sigmas,
                   const vec &extra_parms,
                   const CharacterVector &families,
                   const CharacterVector &links,
                   const field<uvec> &idL_lp_fast,
                   const field<vec> &prior_mean_betas_nHC,
                   const field<mat> &prior_U_Sigma_betas_nHC,
                   const field<mat> &chol_vcov_prop_betas,
                   const field<uvec> &ind_FE_nHC,
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
                   const uword &RMu_it_thld, // robins monro univariate iterations threshold to update the scale. Updates scale from it > RMu_it_thld
                   const double &RMu_tacce, // robins monro univariate target acceptance
                   const field<uvec> &unq_idL,
                   const uword &n_b, // number of unique subjects
                   const uword &n_outcomes
) {
 
 // update eta
 eta = linpred_mixed(X, betas, Z, b, idL);
 
 // update logLik_surv
 Wlong_H = calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, betas, b, 
                               id_H_, FunForms, FunForms_ind); 
 WlongH_alphas = Wlong_H * alphas;

 if (any_event) {
   Wlong_h = calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, betas, b,
                             id_h, FunForms, FunForms_ind);
   Wlongh_alphas = Wlong_h * alphas;
 }
 
 if (any_interval) {
   Wlong_H2 = calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, betas, b, 
                              id_H_, FunForms, FunForms_ind);
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

     uvec ind_j = ind_FE_nHC.at(j);
     // denominator
     double sum_logLik_long_j = sum( log_long_i(y.at(j), eta.at(j), sigmas.at(j), 
                                                extra_parms.at(j), std::string(families[j]), 
                                                std::string(links[j]), idL_lp_fast.at(j))
     );
     
     double logLik_FE_j = log_dmvnrm_chol(betas.at(j).rows(ind_j) - prior_mean_betas_nHC.at(j), 
                                          prior_U_Sigma_betas_nHC.at(j)).at(0); 

     double denominator_j = sum_logLik_long_j + sum(logLik_surv) + logLik_FE_j;
     
     // proposal
     field<vec> betas_prop = betas;
     betas_prop.at(j).rows(ind_j) =+ propose_mvnorm_vec(1, chol_vcov_prop_betas.at(j),
                                                        scale_betas.at(j)); // vec

     // logLik_prior proposal
     double logLik_FE_j_prop = log_dmvnrm_chol(betas_prop.at(j).rows(ind_j) - prior_mean_betas_nHC.at(j), 
                                               prior_U_Sigma_betas_nHC.at(j)).at(0); // vec 1x1
     
     // logLik_long proposal
     field<vec> eta_prop = linpred_mixed(X, betas_prop, Z, b, idL);
     /*
      * To improve: write the function linpred_mixed_i() to work with one outcome
      */ 
     
     double sum_logLik_long_j_prop = sum( log_long_i(y.at(j), eta_prop.at(j), sigmas.at(j), 
                                                     extra_parms.at(j), std::string(families[j]), 
                                                     std::string(links[j]), idL_lp_fast.at(j))
     );
     
     // logLik_surv proposal
     mat Wlong_H_prop = calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, betas_prop, b, 
                                        id_H_, FunForms, FunForms_ind); 
     vec WlongH_alphas_prop = Wlong_H_prop * alphas;
     
     
     mat Wlong_h_prop;
     vec Wlongh_alphas_prop;
     if (any_event) {
       Wlong_h_prop = calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, betas_prop, b,
                                      id_h, FunForms, FunForms_ind);
       Wlongh_alphas_prop = Wlong_h_prop * alphas;
     }
     
     mat Wlong_H2_prop;
     vec WlongH2_alphas_prop;
     if (any_interval) {
       Wlong_H2_prop = calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, betas_prop, b, 
                                       id_H_, FunForms, FunForms_ind);
       WlongH2_alphas_prop = Wlong_H2_prop * alphas;
     }
     
     vec logLik_surv_prop = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas, 
                                     WH_gammas, Wh_gammas, WH2_gammas, 
                                     WlongH_alphas_prop, Wlongh_alphas_prop, WlongH2_alphas_prop, 
                                     log_Pwk, log_Pwk2, id_H_fast, id_h_fast, 
                                     which_event, which_right_event, which_left, 
                                     any_interval, which_interval);
     
     // numerator
     double numerator_j = sum_logLik_long_j_prop + sum(logLik_surv_prop) + logLik_FE_j_prop;
     
     // hastings ratio
     double log_ratio_j = numerator_j - denominator_j;
     
     if(std::isfinite(log_ratio_j) && std::exp(log_ratio_j) > R::runif(0.0, 1.0)) {
       
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
       
       if(it > RMu_it_thld) {
         scale_betas.at(j) = robbins_monro(scale_betas.at(j), acceptance_betas.at(it, j), it, RMu_tacce);
       }
       
     }
     
   }
   
 }
 
 // update logLik_long with all betas
 logLik_long = log_long(y, eta, sigmas, extra_parms,
                        families, links, idL_lp_fast, unq_idL, n_b);
 
 // save all results
 res_betas.row(it) = docall_rbindF(betas).t();

}

// [[Rcpp::export]]
List update_betas_MH (List model_data, List model_info, List initial_values,
                      List priors, List control, List vcov_prop) {

  // model_info
  field<uvec> FunForms = List2Field_uvec(as<List>(model_info["FunForms_cpp"]), true);
  field<uvec> FunForms_ind = List2Field_uvec(as<List>(model_info["FunForms_ind"]), true);
  
  // model_data
  field<mat> X = List2Field_mat(as<List>(model_data["X"]));
  field<mat> Z = List2Field_mat(as<List>(model_data["Z"]));
  field<mat> y = List2Field_mat(as<List>(model_data["y"]));
  field<uvec> idL = List2Field_uvec(as<List>(model_data["idL"]), true);
  vec extra_parms = as<vec>(model_data["extra_parms"]);
  CharacterVector families = as<CharacterVector>(model_info["family_names"]);
  CharacterVector links = as<CharacterVector>(model_info["links"]);
  field<uvec> idL_lp = List2Field_uvec(as<List>(model_data["idL_lp"]), true);
  field<uvec> unq_idL = List2Field_uvec(as<List>(model_data["unq_idL"]), true);
  mat W0_H = as<mat>(model_data["W0_H"]);
  mat W0_h = as<mat>(model_data["W0_h"]);
  mat W0_H2 = as<mat>(model_data["W0_H2"]);
  uvec which_event = as<uvec>(model_data["which_event"]) - 1;
  uvec which_interval = as<uvec>(model_data["which_interval"]) - 1;
  bool any_gammas = as<bool>(model_data["any_gammas"]);
  mat W_H = as<mat>(model_data["W_H"]);
  mat W_h = as<mat>(model_data["W_h"]);
  mat Wlong_H = docall_cbindL(as<List>(model_data["Wlong_H"]));
  mat Wlong_h = docall_cbindL(as<List>(model_data["Wlong_h"]));
  mat Wlong_H2 = docall_cbindL(as<List>(model_data["Wlong_H2"]));
  vec log_Pwk = as<vec>(model_data["log_Pwk"]);
  vec log_Pwk2 = as<vec>(model_data["log_Pwk2"]);
  uvec id_H = as<uvec>(model_data["id_H"]) - 1;
  uvec id_h = as<uvec>(model_data["id_h"]) - 1;
  uvec which_right = as<uvec>(model_data["which_right"]) - 1;
  uvec which_right_event = join_cols(which_event, which_right);
  uvec which_left = as<uvec>(model_data["which_left"]) - 1;
  uvec has_tilde_betas = as<uvec>(model_data["has_tilde_betas"]);
  field<uvec> ind_FE_nHC = List2Field_uvec(as<List>(model_data["ind_FE_nHC"]), true);
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
  uvec id_H_ = as<uvec>(model_data["id_H_"]) - 1;
  
  // control
  uword n_iter = as<uword>(control["n_iter"]);
  
  // initial values
  field<vec> betas = List2Field_vec(as<List>(initial_values["betas"]));
  field<mat> b = List2Field_mat(as<List>(initial_values["b"]));
  vec sigmas = exp(as<vec>(initial_values["log_sigmas"]));
  vec bs_gammas = as<vec>(initial_values["bs_gammas"]);
  vec alphas = as<vec>(initial_values["alphas"]);
  vec gammas = as<vec>(initial_values["gammas"]);
  
  // priors
  field<vec> prior_mean_betas_nHC = List2Field_vec(as<List>(priors["mean_betas_nHC"]));
  
  //vcov_prop
  field<mat> vcov_prop_betas = List2Field_mat(as<List>(vcov_prop["vcov_prop_betas"]));
  field<mat> chol_vcov_prop_betas = vcov_prop_betas;
  for (uword i = 0; i < chol_vcov_prop_betas.n_elem; ++i) {
    chol_vcov_prop_betas.at(i) = chol(vcov_prop_betas.at(i));
  }
  
  // other
  vec betas_vec = docall_rbindF(betas);
  uword n_betas = betas_vec.n_rows;
  mat res_betas(n_iter, n_betas, fill::zeros);
  uword n_outcomes = y.n_elem;
  mat acceptance_betas(n_iter, n_outcomes, fill::zeros);
  vec scale_betas = create_init_scale(n_outcomes);
  
  //
  field<vec> eta = linpred_mixed(X, betas, Z, b, idL);
  mat b_mat = docall_cbindF(b);
  uword n_b = b_mat.n_rows;
  field<uvec> idL_lp_fast(idL_lp.n_elem);
  for (uword i = 0; i < idL_lp.n_elem; ++i) {
    idL_lp_fast.at(i) = create_fast_ind(idL_lp.at(i) + 1);
  }
  
  vec logLik_long = log_long(y, eta, sigmas, extra_parms, families, links,
                             idL_lp_fast, unq_idL, n_b);
  
  // preliminaries
  uvec id_H_fast = create_fast_ind(id_H + 1);
  uvec id_h_fast = create_fast_ind(id_h + 1);
  bool any_event = which_event.n_rows > 0;
  bool any_interval = which_interval.n_rows > 0;
  
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
  vec WlongH2_alphas(W0_H2.n_rows);
  if (any_event) {
    Wlongh_alphas = Wlong_h * alphas;
  }
  if (any_interval) {
    WlongH2_alphas = Wlong_H2 * alphas;
  }
  
  vec logLik_surv = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                             WH_gammas, Wh_gammas, WH2_gammas,
                             WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                             log_Pwk, log_Pwk2, id_H_fast, id_h_fast,
                             which_event, which_right_event, which_left,
                             any_interval, which_interval);

  field<mat> prior_Tau_betas_nHC = List2Field_mat(as<List>(priors["Tau_betas_nHC"]));
  // ?? Add the following to mcmc_fit.cpp
  /*
   * 1. Calculating the chol(vcov) avoids repetead (unnecessary) chol() calculations for each outcome
   * at each iteration inside the update_betas(). 
   * 2. When using log_dmvnrm_chol(), instead of log_dmvnrm(), we only require (up to) 
   * one chol() per outcome.
   * 3. log_dmvnrm_chol() expects a chol(Sigma) and not a chol(Tau).
   * 4. When using inv(diagmat()) I am assuming that prior_Tau_FE.at(j) is diagonal.
   */
  field<mat> prior_U_Sigma_betas_nHC(prior_Tau_betas_nHC.n_elem); 
  for (uword i = 0; i < prior_U_Sigma_betas_nHC.n_elem; ++i) {
    if (!has_tilde_betas.at(i)) {continue;} //?? when adding to mcmc_fit.cpp check if has_tilde_betas is defined before
    prior_U_Sigma_betas_nHC.at(i) = chol( inv(diagmat( prior_Tau_betas_nHC.at(i) )) );
  }

  //////////////////////////////////////////////////////////////////////////////

  for (uword it = 0; it < n_iter; ++it) {

    update_betas (betas,
                  res_betas,
                  acceptance_betas,
                  scale_betas,
                  eta,
                  logLik_long,
                  logLik_surv,
                  Wlong_H, Wlong_h, Wlong_H2,
                  WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                  it, // it
                  has_tilde_betas,
                  X, 
                  Z,
                  b,
                  idL,
                  y,
                  sigmas,
                  extra_parms,
                  families,
                  links,
                  idL_lp_fast,
                  prior_mean_betas_nHC,
                  prior_U_Sigma_betas_nHC,
                  chol_vcov_prop_betas,
                  ind_FE_nHC,
                  X_H, X_h, X_H2,
                  Z_H, Z_h, Z_H2,
                  U_H, U_h, U_H2,
                  Wlong_bar,
                  id_H_, id_h,
                  FunForms,
                  FunForms_ind,
                  alphas,
                  any_event, any_interval,
                  W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                  WH_gammas, Wh_gammas, WH2_gammas,
                  log_Pwk, log_Pwk2,
                  id_H_fast, id_h_fast,
                  which_event, which_right_event, which_left,
                  which_interval,
                  20, //RMu_it_thld, // robins monro univariate iterations threshold to update the scale. Updates scale from it > RMu_it_thld
                  0.234, //RMu_tacce, // robins monro univariate target acceptance
                  unq_idL,
                  n_b, // number of unique subjects
                  n_outcomes);
  
  }
  
  return List::create(
    Named("betas") = betas,
    Named("res_betas") = res_betas,
    Named("acceptance_betas") = acceptance_betas
  );

}
//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////       FUNS /////////////////////////////////////////////
static double const Const_Unif_Proposal = 0.5 * std::pow(12.0, 0.5);
static double const log2pi = std::log(2.0 * M_PI);

double robbins_monro (const double &scale, const double &acceptance_it,
                      const int &it, const double &target_acceptance = 0.45) {
  double step_length = scale / (target_acceptance * (1 - target_acceptance));
  double out;
  if (acceptance_it > 0) {
    out = scale + step_length * (1 - target_acceptance) / (double)it;
  } else {
    out = scale - step_length * target_acceptance / (double)it;
  }
  return out;
}

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

void inplace_LowerTrimat_mult (rowvec &x, const mat &trimat) {
  // in-place multiplication of x with an lower triangular matrix trimat
  // because in-place assignment is much faster but careful in use because
  // it changes the input vector x, i.e., not const
  uword const n = trimat.n_cols;
  for (uword j = 0; j < n; ++j) {
    double tmp(0.0);
    for (uword i = j; i < n; ++i)
      tmp += trimat.at(i, j) * x.at(i);
    x.at(j) = tmp;
  }
}

mat cov2cor (const mat &V) {
  vec Is = sqrt(1 / V.diag());
  mat out = V.each_col() % Is;
  out.each_row() %= Is.t();
  return out;
}

mat cor2cov (const mat &R, const vec &sds) {
  mat out = R.each_col() % sds;
  out.each_row() %= sds.t();
  return out;
}

vec group_sum (const vec &x, const uvec &ind) {
  vec cumsum_x = cumsum(x);
  vec out = cumsum_x.elem(ind);
  out.insert_rows(0, 1);
  out = diff(out);
  return out;
}

vec create_init_scale(const uword &n, const double &fill_val = 0.1) {
  vec out(n);
  out.fill(fill_val);
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

field<mat> mat2field_mat (const mat &b, const field<uvec> &ind_RE) {
  uword n = ind_RE.n_elem;
  field<mat> out(n);
  for (uword i = 0; i < n; i++) {
    out.at(i) = b.cols(ind_RE.at(i));
  }
  return out;
}

field<mat> create_storage(const field<vec> &F, const int &n_iter) {
  uword n = F.size();
  field<mat> out(n);
  for (uword i = 0; i < n; ++i) {
    vec aa = F.at(i);
    uword n_i = aa.n_rows;
    mat tt(n_iter, n_i, fill::zeros);
    out.at(i) = tt;
  }
  return out;
}

vec Wlong_alphas_fun (const field<mat> &Mats, const field<vec> &coefs) {
  uword n = Mats.n_elem;
  uword n_rows = Mats.at(0).n_rows;
  vec out(n_rows, fill::zeros);
  for (uword k = 0; k < n; ++k) {
    out += Mats.at(k) * coefs.at(k);
  }
  return out;
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

uvec create_fast_ind (const uvec &group) {
  uword l = group.n_rows;
  uvec ind = find(group.rows(1, l - 1) != group.rows(0, l - 2));
  uword k = ind.n_rows;
  ind.insert_rows(k, 1);
  ind.at(k) = l - 1;
  return ind;
}

double logPrior(const vec &x, const vec &mean, const mat &Tau,
                const double tau = 1.0) {
  vec z = (x - mean);
  double out = - 0.5 * tau * as_scalar(z.t() * Tau * z);
  return out;
}

vec propose_norm (const vec &thetas, const vec &scale, const uword &i) {
  vec proposed_thetas = thetas;
  proposed_thetas.at(i) = R::rnorm(thetas.at(i), scale.at(i));
  return proposed_thetas;
}

vec propose_unif (const vec &thetas, const vec &scale, const uword &i) {
  vec proposed_thetas = thetas;
  proposed_thetas.at(i) = R::runif(thetas.at(i) - Const_Unif_Proposal * scale.at(i),
                     thetas.at(i) + Const_Unif_Proposal * scale.at(i));
  return proposed_thetas;
}

vec propose_lnorm (const vec &thetas, const double &log_mu_i, const vec &scale,
                   const uword &i) {
  vec proposed_thetas = thetas;
  proposed_thetas.at(i) = R::rlnorm(log_mu_i, scale.at(i));
  return proposed_thetas;
}

vec propose_norm_mala (const vec &thetas, const vec &scale,
                       const double &deriv, const uword &i) {
  vec proposed_thetas = thetas;
  double mu = thetas.at(i) + 0.5 * scale.at(i) * deriv;
  double sigma = sqrt(scale.at(i));
  proposed_thetas.at(i) = R::rnorm(mu, sigma);
  return proposed_thetas;
}

field<vec> propose_field (const field<vec>& thetas,
                          const field<vec>& scale,
                          const uword &k, const uword &i) {
  field<vec> proposed_thetas = thetas;
  proposed_thetas.at(k).at(i) = R::rnorm(thetas.at(k).at(i),
                     scale.at(k).at(i));
  return proposed_thetas;
}

mat rnorm_mat (const uword& rows, const uword& cols) {
  mat out(rows, cols);
  out.each_col([&](vec& x) {x = as<vec>(rnorm(rows)); } );
  return out;
}

// S is the Cholesky factorisation of vcov_prep_RE which needs to be doen outside MCMC loop
// currently with rnorm_mat but we need to check if sth changed with the seeds in Armadillo
// maybe we can go back to randn() [faster]
cube propose_mvnorm_cube (const int& n, const cube& S, const vec& scale) {
  uword ncol_per_slice = S.n_cols;
  uword slices = S.n_slices;
  cube out(n, ncol_per_slice, slices);
  for (uword i = 0; i < slices; i++) {
    out.slice(i) = scale.at(i) * (rnorm_mat(n, ncol_per_slice) * S.slice(i));
  }
  return out;
}

// returns a mat transposed version: same dimensions as b_mat
mat propose_mvnorm_mat (const int& n, const cube& S, const vec& scale) {
  uword ncol_per_slice = S.n_cols;
  uword slices = S.n_slices;
  cube tmp(n, ncol_per_slice, slices);
  for (uword i = 0; i < slices; i++) {
    tmp.slice(i) = scale.at(i) * (rnorm_mat(n, ncol_per_slice) * S.slice(i));
  }
  mat out = tmp.row(0);
  return out.t();
}

vec mu_fun (const vec &eta, const std::string &link) {
  uword n = eta.n_rows;
  vec out(n);
  if (link == "identity") {
    out = eta;
  } else if (link == "logit") {
    out = 1 / (1 + trunc_exp(- eta));
  } else if (link == "probit") {
    out = normcdf(eta);
  } else if (link == "cloglog") {
    out = - trunc_exp(- trunc_exp(eta)) + 1;
  } else if (link == "log") {
    out = trunc_exp(eta);
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

vec log_dbbinom (const vec &x, const vec &size, const vec &prob,
                 const double &phi) {
  vec A = phi * prob;
  vec B = phi * (1 - prob);
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

vec log_dnorm (const vec &x, const vec &mu, const double &sigma) {
  vec sigmas(x.n_rows);
  sigmas.fill(sigma);
  vec out = log_normpdf(x, mu, sigmas);
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

vec log_dgamma (const vec &x, const vec &shape, const vec &scale) {
  uword n = x.n_rows;
  vec out(n);
  for (uword i = 0; i < n; ++i) {
    out.at(i) = R::dgamma(x.at(i), shape.at(i), scale.at(i), 1);
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

vec log_dmvnrm (const mat &x, const mat &D) {
  // fast log density of the multivariate normal distribution
  // D is the covariance matrix.
  using arma::uword;
  uword const n = x.n_rows, xdim = x.n_cols;
  vec out(n);
  mat V = inv(trimatu(chol(D)));
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

field<vec> linpred_mixed_Zb (const field<mat>& Xbetas,
                             const field<mat> &Z, const field<mat> &b,
                             const field<uvec> &id) {
  uword n_outcomes = Z.n_elem;
  field<vec> out(n_outcomes);
  for (uword i = 0; i < n_outcomes; ++i) {
    mat Xbetas_i = Xbetas.at(i);
    mat Z_i = Z.at(i);
    mat b_i = b.at(i);
    uvec id_i = id.at(i);
    out.at(i) = Xbetas_i + arma::sum(Z_i % b_i.rows(id_i), 1);
  }
  return out;
}

arma::field<arma::mat> calculate_u(arma::field<arma::mat> Xhc,
                                   arma::field<arma::uvec> columns_HC,
                                   arma::field<arma::vec> betas,
                                   arma::field<arma::mat> b,
                                   arma::field<arma::uvec> unq_idL) {
  arma::field<arma::mat>u(b);
  uword n = Xhc.n_elem;
  arma::mat Xhc_i;
  arma::uvec columns_HC_i;
  arma::vec betas_i;
  arma::mat b_i;
  arma::uvec unq_idL_i;
  arma::mat mean_b_i;
  uword ncol_b_i;
  arma::uvec index;
  arma::uvec cindex;
  for (uword i = 0; i < n; i++) {
    Xhc_i = Xhc(i);
    columns_HC_i = columns_HC(i);
    betas_i = betas(i);
    b_i = b(i);
    unq_idL_i = unq_idL(i);
    mean_b_i = b_i * 0;
    ncol_b_i = b_i.n_cols;
    for (uword j = 0; j < ncol_b_i; j++) {
      index = find(columns_HC_i == j + 1);
      cindex = j;
      mean_b_i(unq_idL_i - 1, cindex) = Xhc_i.cols(index) * betas_i(index);
    }
    u(i) = b_i + mean_b_i;
  }
  return(u);
}

arma::field<arma::mat> calculate_u_mean(arma::field<arma::mat> Xhc,
                                        arma::field<arma::uvec> columns_HC,
                                        arma::field<arma::vec> betas,
                                        arma::field<arma::mat> b,
                                        arma::field<arma::uvec> unq_idL) {
  arma::field<arma::mat>u(b);
  uword n = Xhc.n_elem;
  arma::mat Xhc_i;
  arma::uvec columns_HC_i;
  arma::vec betas_i;
  arma::mat b_i;
  arma::uvec unq_idL_i;
  arma::mat mean_b_i;
  uword ncol_b_i;
  arma::uvec index;
  arma::uvec cindex;
  for (uword i = 0; i < n; i++) {
    Xhc_i = Xhc(i);
    columns_HC_i = columns_HC(i);
    betas_i = betas(i);
    b_i = b(i);
    unq_idL_i = unq_idL(i);
    mean_b_i = b_i * 0;
    ncol_b_i = b_i.n_cols;
    for (uword j = 0; j < ncol_b_i; j++) {
      index = find(columns_HC_i == j + 1);
      cindex = j;
      mean_b_i(unq_idL_i - 1, cindex) = Xhc_i.cols(index) * betas_i(index);
    }
    u(i) = mean_b_i;
  }
  return(u);
}

field<mat> Xbeta_calc (const field<mat> &X, const field<vec> &betas) {
  uword n = X.n_elem;
  field<mat> out(n);
  for (uword i = 0; i < n; i++) {
    out.at(i) = X.at(i) * betas.at(i);
  }
  return out;
}

cube chol_cube (const cube& S) {
  cube out = S;
  out.each_slice( [&] (mat& X) {chol(X); } );
  return out;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  LOG DENS /////////////////////////////////////////////

vec log_long (const field<mat> &y, const field<vec> &eta, const vec &sigmas,
              const vec &extra_parms, const CharacterVector &families,
              const CharacterVector &links, const field<uvec> &ids,
              const field<uvec> &unq_ids) {
  uword n_outcomes = y.size();
  uvec ns(n_outcomes);
  for (uword i = 0; i < n_outcomes; ++i) {
    ns.at(i) = ids.at(i).n_rows;
  }
  uword n = ns.max();
  vec out(n, fill::zeros);
  for (uword i = 0; i < n_outcomes; ++i) {
    uvec id_i = ids.at(i);
    uvec unq_id_i = unq_ids.at(i);
    mat y_i = y.at(i);
    uword N = y_i.n_rows;
    vec log_contr(N);
    vec mu_i = mu_fun(eta.at(i), std::string(links[i]));
    double sigma_i = sigmas.at(i);
    double extr_prm_i = extra_parms.at(i);
    if (families[i] == "gaussian") {
      log_contr = log_dnorm(y_i, mu_i, sigma_i);
    } else if (families[i] == "Student-t") {
      log_contr = log_dt((y_i - mu_i) / sigma_i, extr_prm_i) - log(sigma_i);
    } else if (families[i] == "beta") {
      log_contr = log_dbeta(y_i, mu_i * sigma_i, sigma_i * (1 - mu_i));
    } else if (families[i] == "Gamma") {
      log_contr = log_dgamma(y_i, square(mu_i) / sigma_i, sigma_i / mu_i);
    } else if (families[i] == "unit Lindley") {
      vec theta = 1 / mu_i - 1;
      vec comp1 = 2 * log(theta) - log(1 + theta);
      vec comp2 = - 3 * log(1 - y_i);
      vec comp3 = - (theta * y_i) / (1 - y_i);
      log_contr = comp1 + comp2 + comp3;
    } else if (families[i] == "binomial") {
      uword k = y_i.n_cols;
      if (k == 2) {
        // y_i.col(1) has been set to y_i.col(0) + y_i.col(1)
        // in jm_fit(), i.e., y_i.col(1) is already the number of trials
        // not the number of failures
        log_contr = log_dbinom(y_i.col(0), y_i.col(1), mu_i);
      } else {
        log_contr = y_i % log(mu_i) + (1 - y_i) % log(1 - mu_i);
      }
    } else if (families[i] == "poisson") {
      log_contr = log_dpois(y_i, mu_i);
    } else if (families[i] == "negative binomial") {
      log_contr = log_dnbinom(y_i, mu_i, sigma_i);
    }  else if (families[i] == "beta binomial") {
      uword k = y_i.n_cols;
      if (k == 2) {
        // y_i.col(1) has been set to y_i.col(0) + y_i.col(1)
        // in jm_fit(), i.e., y_i.col(1) is already the number of trials
        // not the number of failures
        log_contr = log_dbbinom(y_i.col(0), y_i.col(1), mu_i, sigma_i);
      } else {
        vec ones(n, fill::ones);
        log_contr = log_dbbinom(y_i, ones, mu_i, sigma_i);
      }
    }
    out.elem(unq_id_i) += group_sum(log_contr, id_i);
  }
  return out;
}

vec log_surv (const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
              const vec &W0H2_bs_gammas, const vec &WH_gammas,
              const vec &Wh_gammas, const vec &WH2_gammas,
              const vec &WlongH_alphas, const vec &Wlongh_alphas,
              const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
              const uvec &indFast_H, const uvec &which_event,
              const uvec &which_right_event, const uvec &which_left,
              const bool &any_interval, const uvec &which_interval) {
  vec lambda_H = W0H_bs_gammas + WH_gammas + WlongH_alphas;
  vec H = group_sum(exp(log_Pwk + lambda_H), indFast_H);
  uword n = H.n_rows;
  vec lambda_h(n);
  lambda_h.elem(which_event) = W0h_bs_gammas.elem(which_event) +
    Wh_gammas.elem(which_event) + Wlongh_alphas.elem(which_event);
  vec out(n);
  out.elem(which_right_event) = - H.elem(which_right_event);
  out.elem(which_event) += lambda_h.elem(which_event);
  out.elem(which_left) = log1p(- exp(- H.elem(which_left)));
  vec lambda_H2(lambda_H.n_rows);
  vec H2(n);
  if (any_interval) {
    lambda_H2 = W0H2_bs_gammas + WH2_gammas + WlongH2_alphas;
    H2 = group_sum(exp(log_Pwk2 + lambda_H2), indFast_H);
    out.elem(which_interval) = - H.elem(which_interval) +
      log(- expm1(- H2.elem(which_interval)));
  }
  return out;
}

vec log_re (const mat &b, const mat &L, const vec &sds) {
  mat chol_Sigma = L.each_row() % sds.t();
  vec out = log_dmvnrm_chol(b, chol_Sigma);
  return out;
}

vec logLik (const field<mat> &y, const field<vec> &eta, const vec &sigmas,
            const vec &extra_parms, const CharacterVector &families,
            const CharacterVector &links, const field<uvec> &ids,
            const field<uvec> &unq_ids,
            const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
            const vec &W0H2_bs_gammas, const vec &WH_gammas,
            const vec &Wh_gammas, const vec &WH2_gammas,
            const vec &WlongH_alphas, const vec &Wlongh_alphas,
            const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
            const uvec &indFast_H, const uvec &which_event,
            const uvec &which_right_event, const uvec &which_left,
            const bool &any_interval, const uvec &which_interval,
            const mat &b_mat, const mat &L, const vec &sds) {
  vec logLik_long = log_long(y, eta, sigmas, extra_parms, families, links, ids,
                             unq_ids);
  vec logLik_surv = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                             WH_gammas, Wh_gammas, WH2_gammas, WlongH_alphas,
                             Wlongh_alphas, WlongH2_alphas, log_Pwk,
                             log_Pwk2, indFast_H, which_event,
                             which_right_event, which_left, any_interval,
                             which_interval);
  vec logLik_re = log_re(b_mat, L, sds);
  vec out = logLik_long + logLik_surv + logLik_re;
  return out;
}

double logLik_prior (const mat &L, const vec &sds,
                     const double &prior_D_sds_df, const double &prior_D_sds_sigma,
                     const double &prior_D_L_etaLKJ,
                     const vec &bs_gammas, const vec &gammas, const vec &alphas,
                     const vec &prior_mean_bs_gammas, const mat &prior_Tau_bs_gammas,
                     const vec &prior_mean_gammas, const mat &prior_Tau_gammas,
                     const vec &prior_mean_alphas, const mat &prior_Tau_alphas,
                     const double &tau_bs_gammas,
                     double prior_A_tau_bs_gammas, double prior_B_tau_bs_gammas) {
  double out(0.0);
  out += sum(log_dht(sds, prior_D_sds_sigma, prior_D_sds_df));
  uword p = L.n_rows;
  double log_p_L(0.0);
  for (unsigned i = 1; i < p; ++i) {
    log_p_L += (p - i - 1.0 + 2.0 * prior_D_L_etaLKJ - 2.0) * log(L.at(i, i));
  }
  out += log_p_L;
  out += logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas,
                  tau_bs_gammas);
  out += logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0);
  out += logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
  
  return out;
}

vec log_b (const field<mat> &Xbetas, const field<mat> &Z, const field<mat> &b, const mat &b_mat,
           const field<uvec> &id, const field<mat> &y, const vec &scales,
           const vec &extra_parms, const CharacterVector &families,
           const CharacterVector &links, const field<uvec> &ids,
           const field<uvec> &unq_ids, const mat &L,
           const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
           const vec &W0H2_bs_gammas, const vec &WH_gammas,
           const vec &Wh_gammas, const vec &WH2_gammas,
           const vec &WlongH_alphas, const vec &Wlongh_alphas,
           const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
           const uvec &indFast_H, const uvec &which_event,
           const uvec &which_right_event, const uvec &which_left,
           const bool &any_interval, const uvec &which_interval) {
  field<vec> eta = linpred_mixed_Zb(Xbetas, Z, b, id);
  vec log_lik_y = log_long(y, eta, scales, extra_parms, families, links, ids, unq_ids);
  vec log_pb = log_dmvnrm_chol(b_mat, L);
  vec log_lik_surv = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                              WH_gammas, Wh_gammas, WH2_gammas,
                              WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                              log_Pwk, log_Pwk2, indFast_H,
                              which_event, which_right_event, which_left,
                              any_interval, which_interval);
  vec logLik = log_lik_y + log_pb + log_lik_surv;
  return logLik;
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////  UPDATE FUNS /////////////////////////////////////////////

void update_bs_gammas (vec &bs_gammas, const vec &gammas, const vec &alphas,
                       vec &W0H_bs_gammas, vec &W0h_bs_gammas, vec &W0H2_bs_gammas,
                       const vec &WH_gammas, const vec &Wh_gammas, const vec &WH2_gammas,
                       const vec &WlongH_alphas, const vec &Wlongh_alphas, const vec &WlongH2_alphas,
                       const vec &log_Pwk, const vec &log_Pwk2, const uvec &id_H,
                       const uvec &which_event, const uvec &which_right_event,
                       const uvec &which_left, const uvec &which_interval,
                       const bool &any_event, const bool &any_interval,
                       const vec &prior_mean_bs_gammas, const mat &prior_Tau_bs_gammas,
                       const double &tau_bs_gammas,
                       const vec &prior_mean_gammas, const mat &prior_Tau_gammas,
                       const vec &prior_mean_alphas, const mat &prior_Tau_alphas,
                       vec &logLik_surv, double &denominator_surv, const uword &it,
                       /////
                       const mat &W0_H, const mat &W0_h, const mat &W0_H2,
                       vec &scale_bs_gammas, mat &acceptance_bs_gammas,
                       mat &res_bs_gammas) {
  for (uword i = 0; i < bs_gammas.n_rows; ++i) {
    vec proposed_bs_gammas = propose_norm(bs_gammas, scale_bs_gammas, i);
    vec proposed_W0H_bs_gammas = W0_H * proposed_bs_gammas;
    vec proposed_W0h_bs_gammas(W0_h.n_rows);
    vec proposed_W0H2_bs_gammas(W0_H2.n_rows);
    if (any_event) {
      proposed_W0h_bs_gammas = W0_h * proposed_bs_gammas;
    }
    if (any_interval) {
      proposed_W0H2_bs_gammas = W0_H2 * proposed_bs_gammas;
    }
    vec logLik_surv_proposed =
      log_surv(proposed_W0H_bs_gammas, proposed_W0h_bs_gammas, proposed_W0H2_bs_gammas,
               WH_gammas, Wh_gammas, WH2_gammas,
               WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
               log_Pwk, log_Pwk2, id_H,
               which_event, which_right_event, which_left,
               any_interval, which_interval);
    double numerator_surv =
      sum(logLik_surv_proposed) +
      logPrior(proposed_bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
      logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
      logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    double log_ratio = numerator_surv - denominator_surv;
    if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
      bs_gammas = proposed_bs_gammas;
      W0H_bs_gammas = proposed_W0H_bs_gammas;
      if (any_event) {
        W0h_bs_gammas = proposed_W0h_bs_gammas;
      }
      if (any_interval) {
        W0H2_bs_gammas = proposed_W0H2_bs_gammas;
      }
      logLik_surv = logLik_surv_proposed;
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

void update_gammas (const vec &bs_gammas, vec &gammas, const vec &alphas,
                    const vec &W0H_bs_gammas, const vec &W0h_bs_gammas, const vec &W0H2_bs_gammas,
                    vec &WH_gammas, vec &Wh_gammas, vec &WH2_gammas,
                    const vec &WlongH_alphas, const vec &Wlongh_alphas, const vec &WlongH2_alphas,
                    const vec &log_Pwk, const vec &log_Pwk2, const uvec &id_H,
                    const uvec &which_event, const uvec &which_right_event,
                    const uvec &which_left, const uvec &which_interval,
                    const bool &any_event, const bool &any_interval,
                    const vec &prior_mean_bs_gammas, const mat &prior_Tau_bs_gammas,
                    const double &tau_bs_gammas,
                    const vec &prior_mean_gammas, const mat &prior_Tau_gammas,
                    const vec &prior_mean_alphas, const mat &prior_Tau_alphas,
                    vec &logLik_surv, double &denominator_surv, const uword &it,
                    /////
                    const mat &W_H, const mat &W_h, const mat &W_H2,
                    vec &scale_gammas, mat &acceptance_gammas, mat &res_gammas) {
  for (uword i = 0; i < gammas.n_rows; ++i) {
    vec proposed_gammas = propose_norm(gammas, scale_gammas, i);
    vec proposed_WH_gammas = W_H * proposed_gammas;
    vec proposed_Wh_gammas(W_h.n_rows);
    vec proposed_WH2_gammas(W_H2.n_rows);
    if (any_event) {
      proposed_Wh_gammas = W_h * proposed_gammas;
    }
    if (any_interval) {
      proposed_WH2_gammas = W_H2 * proposed_gammas;
    }
    vec logLik_surv_proposed =
      log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
               proposed_WH_gammas, proposed_Wh_gammas, proposed_WH2_gammas,
               WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
               log_Pwk, log_Pwk2, id_H,
               which_event, which_right_event, which_left,
               any_interval, which_interval);
    double numerator_surv =
      sum(logLik_surv_proposed) +
      logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
      logPrior(proposed_gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
      logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    double log_ratio = numerator_surv - denominator_surv;
    if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
      gammas = proposed_gammas;
      WH_gammas = proposed_WH_gammas;
      if (any_event) {
        Wh_gammas = proposed_Wh_gammas;
      }
      if (any_interval) {
        WH2_gammas = proposed_WH2_gammas;
      }
      logLik_surv = logLik_surv_proposed;
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

void update_alphas (const vec &bs_gammas, const vec &gammas, vec &alphas,
                    const vec &W0H_bs_gammas, const vec &W0h_bs_gammas, const vec &W0H2_bs_gammas,
                    const vec &WH_gammas, const vec &Wh_gammas, const vec &WH2_gammas,
                    vec &WlongH_alphas, vec &Wlongh_alphas, vec &WlongH2_alphas,
                    const vec &log_Pwk, const vec &log_Pwk2, const uvec &id_H,
                    const uvec &which_event, const uvec &which_right_event,
                    const uvec &which_left, const uvec &which_interval,
                    const bool &any_event, const bool &any_interval,
                    const vec &prior_mean_bs_gammas, const mat &prior_Tau_bs_gammas,
                    const double &tau_bs_gammas,
                    const vec &prior_mean_gammas, const mat &prior_Tau_gammas,
                    const vec &prior_mean_alphas, const mat &prior_Tau_alphas,
                    vec &logLik_surv, double &denominator_surv, const uword &it,
                    /////
                    const mat &Wlong_H, const mat &Wlong_h, const mat &Wlong_H2,
                    vec &scale_alphas, mat &acceptance_alphas, mat &res_alphas) {
  for (uword i = 0; i < alphas.n_rows; ++i) {
    vec proposed_alphas = propose_norm(alphas, scale_alphas, i);
    vec proposed_WlongH_alphas = Wlong_H * proposed_alphas;
    vec proposed_Wlongh_alphas(Wlong_h.n_rows);
    if (any_event) {
      proposed_Wlongh_alphas = Wlong_h * proposed_alphas;
    }
    vec proposed_WlongH2_alphas(Wlong_H2.n_rows);
    if (any_interval) {
      proposed_WlongH2_alphas = Wlong_H2 * proposed_alphas;
    }
    vec logLik_surv_proposed =
      log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
               WH_gammas, Wh_gammas, WH2_gammas,
               proposed_WlongH_alphas, proposed_Wlongh_alphas, proposed_WlongH2_alphas,
               log_Pwk, log_Pwk2, id_H,
               which_event, which_right_event, which_left,
               any_interval, which_interval);
    double numerator_surv =
      sum(logLik_surv_proposed) +
      logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
      logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
      logPrior(proposed_alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    double log_ratio = numerator_surv - denominator_surv;
    if (std::isfinite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
      alphas = proposed_alphas;
      WlongH_alphas = proposed_WlongH_alphas;
      if (any_event) {
        Wlongh_alphas = proposed_Wlongh_alphas;
      }
      if (any_interval) {
        WlongH2_alphas = proposed_WlongH2_alphas;
      }
      logLik_surv = logLik_surv_proposed;
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

void update_Wlong (mat &Wlong_H, mat &Wlong_h, mat &Wlong_H2,
                   const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
                   const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
                   const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
                   const field<vec> &betas, const field<mat> &b,
                   const uvec &id_H, const uvec &id_h,
                   const field<uvec> &FunForms, const field<uvec> &FunForms_ind,
                   const bool &any_event, const bool &any_interval) {
  field<mat> eta_H = linpred_surv(X_H, betas, Z_H, b, id_H);
  Wlong_H = docall_cbindF(create_Wlong(eta_H, FunForms, U_H, FunForms_ind));
  if (any_event) {
    field<mat> eta_h = linpred_surv(X_h, betas, Z_h, b, id_h);
    Wlong_h = docall_cbindF(create_Wlong(eta_h, FunForms, U_h, FunForms_ind));
  }
  if (any_interval) {
    field<mat> eta_H2 = linpred_surv(X_H2, betas, Z_H2, b, id_H);
    Wlong_H2 = docall_cbindF(create_Wlong(eta_H2, FunForms, U_H2, FunForms_ind));
  }
}

void update_mean_u (field<mat> &mean_u, const field<vec> &betas,
                    const field<mat> &Xbase, const field<uvec> &x_in_z,
                    const field<uvec> &baseline, const field<uvec> &unq_idL) {
  uword n = mean_u.n_elem;
  for (uword i = 0; i < n; ++i) {
    vec betas_i = betas.at(i);
    mat Xbase_i = Xbase.at(i);
    uvec xinz_i = x_in_z.at(i);
    uvec base_i = baseline.at(i);
    uvec rowind_i = unq_idL.at(i);
    uword k = xinz_i.n_rows;
    if (mean_u.at(i).n_cols == k) {
      mean_u.at(i).each_row() = betas_i.rows(xinz_i).t();
    } else {
      mean_u.at(i).cols(0, k - 1).each_row() = betas_i.rows(xinz_i).t();
    }
    if (is_finite(Xbase_i)) {
      mean_u.at(i)(rowind_i) = Xbase_i * betas_i.rows(base_i);
    }
  }
}

vec log_dht (const vec &x, const double &sigma = 10.0,
             const double &df = 3.0) {
  // log density of half Student's t with scale sigma and df degrees of freedom
  // https://en.m.wikipedia.org/wiki/Folded-t_and_half-t_distributions
  uword n = x.n_rows;
  vec out(n);
  double log_const = std::log(2.0) + lgamma(0.5 * (df + 1)) - lgamma(0.5 * df) -
    0.5 * (std::log(df) + std::log(M_PI)) - std::log(sigma);
  vec log_kernel = - 0.5 * (df + 1.0) * log(1.0 + square(x) / (df * pow(sigma, 2.0)));
  out = log_const + log_kernel;
  return out;
}

double logPC_D_sds (const vec &sds, const mat &L, const mat &b,
                    const double &prior_D_sds_df,
                    const double &prior_D_sds_sigma) {
  mat chol_Sigma = L.each_row() % sds.t();
  double log_p_b = sum(log_dmvnrm_chol(b, chol_Sigma));
  double log_p_sds = sum(log_dht(sds, prior_D_sds_sigma, prior_D_sds_df));
  double out = log_p_b + log_p_sds;
  return out;
}

double logPC_D_L (const mat &L, const vec &sds, const mat &b,
                  const double &prior_D_L_etaLKJ) {
  uword p = L.n_rows;
  mat chol_Sigma = L.each_row() % sds.t();
  double log_p_b = sum(log_dmvnrm_chol(b, chol_Sigma));
  double log_p_L(0.0);
  for (unsigned i = 1; i < p; ++i) {
    log_p_L += (p - i - 1.0 + 2.0 * prior_D_L_etaLKJ - 2.0) * log(L.at(i, i));
  }
  double out = log_p_b + log_p_L;
  return out;
}

double deriv_L (const mat &L, const vec &sds, const mat &b,
                const double &log_target, const uword &i,
                const uvec &upper_part,
                const double &prior_D_L_etaLKJ,
                const char &direction = 'b', const double &eps = 1e-06) {
  uword n = L.n_rows;
  uword upper_part_i = upper_part.at(i);
  uword column = floor(upper_part_i / n);
  mat L_eps = L;
  if (direction == 'f') {
    L_eps(upper_part_i) += L_eps(upper_part_i) * eps;
  } else {
    L_eps(upper_part_i) -= L_eps(upper_part_i) * eps;
  }
  vec ll = L_eps.submat(0, column, column - 1, column);
  double ss = dot(ll, ll);
  if (ss > 1) return datum::nan;
  L_eps.at(column, column) = sqrt(1 - ss);
  double out(0.0);
  if (direction == 'f') {
    out = (logPC_D_L(L_eps, sds, b, prior_D_L_etaLKJ) - log_target) / eps;
  } else {
    out = (log_target - logPC_D_L(L_eps, sds, b, prior_D_L_etaLKJ)) / eps;
  }
  return out;
}

mat propose_L (const mat &L, const vec &scale, const uvec &upper_part,
               const double &deriv, const uword &i, const bool &mala = false) {
  mat proposed_L = L;
  vec l = L(upper_part);
  vec proposed_l(l.n_rows);
  if (mala) {
    if (std::isfinite(deriv)) {
      proposed_l = propose_norm_mala(l, scale, deriv, i);
    } else {
      return proposed_L.fill(datum::nan);
    }
  } else {
    proposed_l = propose_unif(l, scale, i);
  }
  proposed_L(upper_part) = proposed_l;
  uword n = L.n_rows;
  uword column = floor(upper_part.at(i) / n);
  vec ll = proposed_L.submat(0, column, column - 1, column);
  double ss = dot(ll, ll);
  if (ss > 1) return proposed_L.fill(datum::nan);
  proposed_L.at(column, column) = sqrt(1 - ss);
  return proposed_L;
}

void update_D (mat &L, vec &sds, const mat &b,
               const uvec &upper_part,
               const double &prior_D_sds_df,
               const double &prior_D_sds_sigma,
               const double &prior_D_L_etaLKJ,
               const int &it, const bool &MALA,
               mat &res_sds, mat &res_L,
               vec &scale_sds, vec &scale_L,
               mat &acceptance_sds, mat &acceptance_L) {
  uword n_sds = sds.n_rows;
  uword n_L = upper_part.n_rows;
  double denominator_sds = logPC_D_sds(sds, L, b, prior_D_sds_df,
                                       prior_D_sds_sigma);
  for (uword i = 0; i < n_sds; ++i) {
    double SS = 0.5 * pow(scale_sds.at(i), 2.0);
    double log_mu_current = log(sds.at(i)) - SS;
    vec proposed_sds = propose_lnorm(sds, log_mu_current, scale_sds, i);
    double numerator_sds = logPC_D_sds(proposed_sds, L, b,
                                       prior_D_sds_df, prior_D_sds_sigma);
    double log_mu_proposed = log(proposed_sds.at(i)) - SS;
    double log_ratio_sds = numerator_sds - denominator_sds +
      R::dlnorm(sds.at(i), log_mu_proposed, scale_sds.at(i), true) -
      R::dlnorm(proposed_sds.at(i), log_mu_current, scale_sds.at(i), true);
    if (std::isfinite(log_ratio_sds) && exp(log_ratio_sds) > R::runif(0.0, 1.0)) {
      sds = proposed_sds;
      denominator_sds = numerator_sds;
      acceptance_sds.at(it, i) = 1;
    }
    if (it > 19) {
      scale_sds.at(i) =
        robbins_monro(scale_sds.at(i), acceptance_sds.at(it, i),
                      it);
    }
    res_sds.at(it, i) = sds.at(i);
  }
  double denominator_L = logPC_D_L(L, sds, b, prior_D_L_etaLKJ);
  for (uword i = 0; i < n_L; ++i) {
    uword upper_part_i = upper_part.at(i);
    double deriv_current(0.0);
    double mu_current(0.0);
    mat proposed_L = L;
    if (MALA) {
      deriv_current = deriv_L(L, sds, b, denominator_L, i, upper_part,
                              prior_D_L_etaLKJ);
      mu_current = L.at(upper_part_i) + 0.5 * scale_L.at(i) * deriv_current;
      proposed_L = propose_L(L, scale_L, upper_part, deriv_current, i, true);
    } else {
      proposed_L = propose_L(L, scale_L, upper_part, deriv_current, i);
    }
    double numerator_L(0.0);
    double deriv_proposed(0.0);
    double mu_proposed(0.0);
    double log_ratio_L(0.0);
    bool finite_L = proposed_L.is_finite();
    if (finite_L) {
      numerator_L = logPC_D_L(proposed_L, sds, b, prior_D_L_etaLKJ);
      if (MALA) {
        deriv_proposed = deriv_L(proposed_L, sds, b, numerator_L,
                                 i, upper_part, prior_D_L_etaLKJ);
        mu_proposed = proposed_L.at(upper_part_i) +
          0.5 * scale_L.at(i) * deriv_proposed;
        log_ratio_L = numerator_L - denominator_L +
          log_normpdf(L.at(upper_part_i), mu_proposed, sqrt(scale_L.at(i))) -
          log_normpdf(proposed_L.at(upper_part_i), mu_current, sqrt(scale_L.at(i)));
      } else {
        log_ratio_L = numerator_L - denominator_L;
      }
    }
    if (finite_L && std::isfinite(log_ratio_L) &&
        exp(log_ratio_L) > R::runif(0.0, 1.0)) {
      L = proposed_L;
      denominator_L = numerator_L;
      acceptance_L.at(it, i) = 1;
    }
    if (it > 19) {
      scale_L.at(i) =
        robbins_monro(scale_L.at(i), acceptance_L.at(it, i),
                      it);
    }
    res_L.at(it, i) = L.at(upper_part_i);
  }
}


void update_b (field<mat> &b, mat &b_mat, const field<mat> &Xbetas, const field<mat> &Z, const field<uvec> &id, 
               const field<mat> &y, const vec &extra_parms, 
               const CharacterVector &families, const CharacterVector &links, const field<uvec> &ids,
               const field<uvec> &unq_ids, const mat &L, 
               const vec &W0H_bs_gammas, const vec &W0h_bs_gammas,
               const vec &W0H2_bs_gammas, const vec &WH_gammas,
               const vec &Wh_gammas, const vec &WH2_gammas,
               const vec &WlongH_alphas, const vec &Wlongh_alphas,
               const vec &WlongH2_alphas, const vec &log_Pwk, const vec &log_Pwk2,
               const uvec &indFast_H, const uvec &which_event,
               const uvec &which_right_event, const uvec &which_left,
               const bool &any_interval, const uvec &which_interval, 
               const field<uvec> ind_RE, const uword &it,
               cube &res_b, vec &scale_b, mat &acceptance_b) {
  mat proposed_b = propose_mvnorm_mat(1, L, scale_b) + b_mat;
  field<mat> proposed_b_field = mat2field_mat(proposed_b, ind_RE); 
  vec numerator_b = log_b(Xbetas, Z, proposed_b_field, proposed_b,
                          id, y, scale_b,
                          extra_parms, families,
                          links, ids,
                          unq_ids, L, 
                          W0H_bs_gammas, W0h_bs_gammas,
                          W0H2_bs_gammas, WH_gammas,
                          Wh_gammas, WH2_gammas,
                          WlongH_alphas, Wlongh_alphas,
                          WlongH2_alphas, log_Pwk, log_Pwk2,
                          indFast_H, which_event,
                          which_right_event, which_left,
                          any_interval, which_interval); 
  vec denominator_b = log_b(Xbetas, Z, b, b_mat,
                            id, y, scale_b,
                            extra_parms, families,
                            links, ids,
                            unq_ids, L, 
                            W0H_bs_gammas, W0h_bs_gammas,
                            W0H2_bs_gammas, WH_gammas,
                            Wh_gammas, WH2_gammas,
                            WlongH_alphas, Wlongh_alphas,
                            WlongH2_alphas, log_Pwk, log_Pwk2,
                            indFast_H, which_event,
                            which_right_event, which_left,
                            any_interval, which_interval); 
  vec log_ratio = numerator_b - denominator_b;
  uword n = log_ratio.n_elem;
  for (uword i = 0; i < n; i++) {
    if (std::isfinite(log_ratio.at(i)) && exp(log_ratio.at(i)) > R::runif(0, 1)) {
      acceptance_b.at(i, it) = 1;
      b_mat.row(i) = proposed_b.row(i);
    }
    if (it > 19) {
      scale_b.at(i) =
        robbins_monro(scale_b.at(i),
                      acceptance_b.at(i, it), it);
    }
    res_b.row(i) = b_mat.row(i);
  }
  b = mat2field_mat(b_mat, ind_RE);
}


// [[Rcpp::export]]
List mcmc_b_test (List model_data, List model_info, List initial_values, 
                  List priors, List control, List vcov_prop) {
  //////////////////////////////////////////////
  field<mat> X = List2Field_mat(as<List>(model_data["X"]));
  field<vec> betas = List2Field_vec(as<List>(initial_values["betas"]));
  field<mat> Z = List2Field_mat(as<List>(model_data["Z"]));
  field<mat> b = List2Field_mat(as<List>(initial_values["b"]));
  mat b_mat = docall_cbindF(b);
  field<vec> y = List2Field_vec(as<List>(initial_values["y"]));
  vec extra_parms = as<vec>(model_data["extra_parms"]);
  CharacterVector families = as<CharacterVector>(model_info["families"]);
  CharacterVector links = as<CharacterVector>(model_info["links"]);
  field<uvec> idL_lp = List2Field_uvec(as<List>(model_data["idL_lp"]));
  field<uvec> unq_idL = List2Field_uvec(as<List>(model_data["unq_idL"]));
  cube S = as<cube>(vcov_prop["vcov_prop_RE"]);
  cube chol_S = chol_cube(S);
  List ind_RE_ = as<List>(model_data["ind_RE"]);
  field<uvec> ind_RE = List2Field_uvec(ind_RE_, true);
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
  vec WlongH2_alphas(W0_H2.n_rows);
  if (any_event) {
    Wlongh_alphas = Wlong_h * alphas;
  }
  if (any_interval) {
    WlongH2_alphas = Wlong_H2 * alphas;
  }
  vec logLik_surv =
    log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
             WH_gammas, Wh_gammas, WH2_gammas,
             WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
             log_Pwk, log_Pwk2, id_H_fast,
             which_event, which_right_event, which_left,
             any_interval, which_interval);
  double denominator_surv =
    sum(logLik_surv) +
    logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
    logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
    logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
  for (uword it = 0; it < n_iter; ++it) {
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
                     logLik_surv, denominator_surv, it,
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
                    logLik_surv, denominator_surv, it,
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
                  logLik_surv, denominator_surv, it,
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
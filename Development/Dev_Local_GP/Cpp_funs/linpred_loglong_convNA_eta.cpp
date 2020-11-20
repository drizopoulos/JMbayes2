//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

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
  } else if (link == "inverse") {
    out = 1 / eta;
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

cube chol_cube (const cube &S) {
  cube out = S;
  out.each_slice([](mat &X){X = chol(X);});
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
    log_contr = log_dbeta(y_i, mu_i * sigma_i, sigma_i * (1 - mu_i));
  } else if (fam_i == "Gamma") {
    log_contr = log_dgamma(y_i, sigma_i, mu_i / sigma_i);
  } else if (fam_i == "unit Lindley") {
    vec theta = 1 / mu_i - 1;
    vec comp1 = 2 * log(theta) - log(1 + theta);
    vec comp2 = - 3 * log(1 - y_i);
    vec comp3 = - (theta * y_i) / (1 - y_i);
    log_contr = comp1 + comp2 + comp3;
  } else if (fam_i == "binomial") {
    uword k = y_i.n_cols;
    if (k == 2) {
      // y_i.col(1) has been set to y_i.col(0) + y_i.col(1)
      // in jm_fit(), i.e., y_i.col(1) is already the number of trials
      // not the number of failures
      log_contr = log_dbinom(y_i.col(0), y_i.col(1), mu_i);
    } else {
      log_contr = y_i % log(mu_i) + (1 - y_i) % log(1 - mu_i);
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

vec log_long (const field<mat> &y, const field<vec> &eta, const vec &sigmas,
              const vec &extra_parms, const CharacterVector &families,
              const CharacterVector &links, const field<uvec> &idFast,
              const field<uvec> &unq_ids) {
  uword n_outcomes = y.size();
  uvec ns(n_outcomes);
  for (uword i = 0; i < n_outcomes; ++i) {
    ns.at(i) = idFast.at(i).n_rows;
  }
  uword n = ns.max();
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
    out.elem(unq_id_i) += log_contr_i;
  }
  return out;
}

// [[Rcpp::export]]
field<vec> fooeta (List model_data, List model_info, List initial_values) {
  field<mat> X = List2Field_mat(as<List>(model_data["X"]));
  field<mat> Z = List2Field_mat(as<List>(model_data["Z"]));
  field<mat> y = List2Field_mat(as<List>(model_data["y"]));
  field<uvec> idL = List2Field_uvec(as<List>(model_data["idL"]), true);
  field<uvec> unq_idL = List2Field_uvec(as<List>(model_data["unq_idL"]), true);
  field<uvec> idL_lp = List2Field_uvec(as<List>(model_data["idL_lp"]), true);
  vec extra_parms = as<vec>(model_data["extra_parms"]);
  field<uvec> ind_RE = List2Field_uvec(as<List>(model_data["ind_RE"]), true);
  CharacterVector families = as<CharacterVector>(model_info["family_names"]);
  CharacterVector links = as<CharacterVector>(model_info["links"]);
  field<vec> betas = List2Field_vec(as<List>(initial_values["betas"]));
  vec sigmas = exp(as<vec>(initial_values["log_sigmas"]));
  field<mat> b = List2Field_mat(as<List>(initial_values["b"]));
  mat b_mat = docall_cbindF(b);
  field<uvec> idL_lp_fast(idL_lp.n_elem);
  for (uword i = 0; i < idL_lp.n_elem; ++i) {
    idL_lp_fast.at(i) = create_fast_ind(idL_lp.at(i));
  }
  field<vec> eta_proposed = linpred_mixed(X, betas, Z, b, idL_lp);
  return eta_proposed;
}
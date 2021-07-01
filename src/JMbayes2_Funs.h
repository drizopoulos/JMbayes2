#ifndef JMBAYES2FUNS_H
#define JMBAYES2FUNS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

static double const Const_Unif_Proposal = 0.5 * std::pow(12.0, 0.5);
static double const log2pi = std::log(2.0 * M_PI);

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
  vec Is = sqrt(1.0 / V.diag());
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
  vec out = cumsum_x.rows(ind);
  out.insert_rows(0, 1);
  out = diff(out);
  return out;
}

vec create_init_scale(const uword &n, const double &fill_val = 0.1) {
  vec out(n);
  out.fill(fill_val);
  return out;
}

field<vec> create_init_scaleF(const field<uvec> &x, const double &fill_val = 0.1) {
  uword n = x.size();
  field<vec> out(n);
  for (uword i = 0; i < n; ++i) {
    uvec x_i = x.at(i);
    uword n_i = x_i.n_rows;
    vec oo(n_i);
    oo.fill(fill_val);
    out.at(i) = oo;
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

field<mat> mat2field (const mat &b, const field<uvec> &ind_RE) {
  uword n = ind_RE.n_elem;
  field<mat> out(n);
  for (uword i = 0; i < n; i++) {
    out.at(i) = b.cols(ind_RE.at(i));
  }
  return out;
}

field<vec> vec2field (const vec &betas, const field<uvec> &ind_FE) {
  uword n = ind_FE.n_elem;
  field<vec> out(n);
  for (uword i = 0; i < n; i++) {
    out.at(i) = betas.rows(ind_FE.at(i));
  }
  return out;
}

field<vec> create_storage (const field<uvec> &F) {
  uword n = F.size();
  field<vec> out(n);
  for (uword i = 0; i < n; ++i) {
    vec tt(F.at(i).n_rows, fill::zeros);
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
  if (l == 1) return group - 1;
  uvec ind = find(group.rows(1, l - 1) != group.rows(0, l - 2));
  uword k = ind.n_rows;
  ind.insert_rows(k, 1);
  ind.at(k) = l - 1;
  return ind;
}

double logPrior(const vec &x, const vec &mean, mat &Tau, const vec &lambda,
                const double &tau, const bool &shrink) {
  vec z = x - mean;
  if (shrink) {
    Tau.diag() = lambda;
  }
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

mat rnorm_mat (const uword &rows, const uword &cols) {
  mat out(rows, cols);
  out.each_col([&](vec& x) {x = as<vec>(rnorm(rows)); } );
  return out;
}

// S is the Cholesky factorisation of vcov_prep_RE which needs to be doen outside MCMC loop
// currently with rnorm_mat but we need to check if sth changed with the seeds in Armadillo
// maybe we can go back to randn() [faster]
cube propose_mvnorm_cube (const int &n, const cube &S, const vec &scale) {
  uword ncol_per_slice = S.n_cols;
  uword slices = S.n_slices;
  cube out(n, ncol_per_slice, slices);
  for (uword i = 0; i < slices; i++) {
    out.slice(i) = scale.at(i) * (rnorm_mat(n, ncol_per_slice) * S.slice(i));
  }
  return out;
}

mat propose_rnorm_mat (const mat &thetas, const mat &scale, const uword &i) {
  mat proposed_thetas = thetas;
  proposed_thetas.col(i) = scale.col(i) % randn(thetas.n_rows, 1) + thetas.col(i);
  return proposed_thetas;
}

mat propose_rnorm_mat2 (const mat &thetas, const mat &scale, const uword &i) {
  mat proposed_thetas = thetas;
  vec out(thetas.n_rows);
  for (uword i = 0; i < thetas.n_rows; i++) {
    out.at(i) = R::rnorm(0.0, 1.0);
  }
  proposed_thetas.col(i) = scale.col(i) % out + thetas.col(i);
  return proposed_thetas;
}

// returns a mat transposed version: same dimensions as b_mat
mat propose_mvnorm_mat (const int &n, const cube &S, const vec &scale) {
  uword ncol_per_slice = S.n_cols;
  uword slices = S.n_slices;
  cube tmp(n, ncol_per_slice, slices);
  for (uword i = 0; i < slices; i++) {
    tmp.slice(i) = scale.at(i) * (rnorm_mat(n, ncol_per_slice) * S.slice(i));
  }
  mat out = tmp.row(0);
  return out.t();
}

vec propose_mvnorm_vec (const mat &U, const double &scale) {
  uword ncols = U.n_cols;
  vec res = scale * trans(rnorm_mat(1, ncols) * U);
  return res;
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

vec log_dht (const vec &x, const vec &sigma, const double &df = 3.0) {
  // log density of half Student's t with scale sigma and df degrees of freedom
  // https://en.m.wikipedia.org/wiki/Folded-t_and_half-t_distributions
  uword n = x.n_rows;
  vec out(n);
  vec log_const = std::log(2.0) + lgamma(0.5 * (df + 1.0)) - lgamma(0.5 * df) -
    0.5 * (std::log(df) + std::log(M_PI)) - log(sigma);
  vec log_kernel = - 0.5 * (df + 1.0) * log(1.0 + square(x) / (df * square(sigma)));
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

field<vec> linpred_mixed_i (const field<vec> eta, const field<mat> &X,
                            const field<vec> &betas, const field<mat> &Z,
                            const field<mat> &b, const field<uvec> &id,
                            const uword &i) {
  field<vec> out = eta;
  out.at(i) = X.at(i) * betas.at(i) +
    arma::sum(Z.at(i) % b.at(i).rows(id.at(i)), 1);
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

void transf_eta (mat &eta, const CharacterVector &fun_nams) {
  uword n = eta.n_cols;
  for (uword i = 0; i < n; i++) {
    if (fun_nams[i] == "identity") continue;
    if (fun_nams[i] == "expit") {
      eta.col(i) = 1.0 / (1.0 + trunc_exp(- eta.col(i)));
    } else if (fun_nams[i] == "exp" || fun_nams[i] == "dexp") {
      eta.col(i) = trunc_exp(eta.col(i));
    } else if (fun_nams[i] == "dexpit") {
      mat pp = 1.0 / (1.0 + trunc_exp(- eta.col(i)));
      eta.col(i) = pp * (1.0 - pp);
    } else if (fun_nams[i] == "log") {
      eta.col(i) = trunc_log(eta.col(i));
    } else if (fun_nams[i] == "log2") {
      eta.col(i) = log2(eta.col(i));
    } else if (fun_nams[i] == "log10") {
      eta.col(i) = log10(eta.col(i));
    } else if (fun_nams[i] == "sqrt") {
      eta.col(i) = sqrt(eta.col(i));
    } else if (fun_nams[i] == "square") {
      eta.col(i) = square(eta.col(i));
    }
  }
}

field<mat> create_Wlong (const field<mat> &eta, const field<uvec> &FunForms,
                         const field<mat> &U, const field<uvec> &ind,
                         const List &Funs_FunForms) {
  uword n_outcomes = eta.n_elem;
  field<mat> out(n_outcomes);
  for (uword i = 0; i < n_outcomes; ++i) {
    CharacterVector Funs_i = Funs_FunForms[i];
    mat eta_i = eta.at(i);
    transf_eta(eta_i, Funs_i);
    uvec FF_i = FunForms.at(i);
    mat U_i = U.at(i);
    uvec ind_i = ind.at(i);
    mat Wlong_i(eta_i.n_rows, U_i.n_cols, fill::ones);
    Wlong_i.cols(FF_i) %= eta_i.cols(ind_i);
    out.at(i) = U_i % Wlong_i;
  }
  return out;
}

mat calculate_Wlong (const field<mat> &X, const field<mat> &Z,
                     const field<mat> &U, const mat &Wlong_bar,
                     const mat &Wlong_sds,
                     const field<vec> &betas, const field<mat> &b,
                     const uvec &id, const field<uvec> &FunForms,
                     const field<uvec> &FunForms_ind,
                     const List &Funs_FunForms) {
  field<mat> eta = linpred_surv(X, betas, Z, b, id);
  mat Wlong =
    docall_cbindF(create_Wlong(eta, FunForms, U, FunForms_ind, Funs_FunForms));
  Wlong.each_row() -= Wlong_bar;
  Wlong.each_row() /= Wlong_sds;
  return Wlong;
}

mat bdiagF (const field<mat> &F) { // builds a block diagonal matrix given a field of matrices
  uword n; n = F.n_elem; // assumes all matrices being square (nrow=ncol), but with different dim
  uword nrows = 0;
  uvec rows(n);
  for (uword i = 0; i < n; i++) {
    rows.at(i) = F.at(i).n_rows;
    nrows += rows.at(i);
  }
  mat B(nrows, nrows, fill::zeros);
  uword ii = 0;
  for (uword i = 0; i < n; i++) {
    B.submat(ii, ii, ii - 1 + rows.at(i), ii - 1 + rows.at(i)) = F.at(i);
    ii += rows.at(i);
  }
  return B;
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

uword n_field (const field<vec> &x) {
  uword n = x.n_rows;
  uword out = 0;
  for (uword i = 0; i < n; ++i)
    out += x.at(i).n_rows;
  return out;
}

field<vec> create_sigmas_field (const field<vec> &sigmas,
                                const uvec &ss_sigmas,
                                const field<uvec> &idL) {
  uword n = sigmas.size();
  field<vec> out(n);
  for (uword i = 0; i < n; ++i) {
    vec sigmas_i = sigmas.at(i);
    uvec id_i = idL.at(i);
    if (ss_sigmas.at(i)) {
      out.at(i) = sigmas_i.rows(id_i);
    } else {
      vec xx(id_i.n_rows);
      xx.fill(as_scalar(sigmas_i));
      out.at(i) = xx;
    }
  }
  return out;
}

vec scalar2vec (const double &x) {
  vec v(1);
  v.fill(x);
  return v;
}

#endif

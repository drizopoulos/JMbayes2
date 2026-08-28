#ifndef JMBAYES2FUNS_H
#define JMBAYES2FUNS_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

static double const Const_Unif_Proposal = 0.5 * std::pow(12.0, 0.5);
static double const log2pi = std::log(2.0 * M_PI);
static double const half_log2pi = 0.5 * std::log(2.0 * M_PI);

double robbins_monro (double scale, double acceptance_it,
                          int it, double target_acceptance = 0.45) {
    double it_d = static_cast<double>(it);
    if (acceptance_it > 0.0) {
        return scale * (1.0 + 1.0 / (target_acceptance * it_d));
    } else {
        return scale * (1.0 - 1.0 / ((1.0 - target_acceptance) * it_d));
    }
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

arma::vec group_sum(const arma::vec& x, const arma::uvec& ind) {
    arma::uword m = ind.n_elem;
    arma::vec out(m);
    // Extract raw pointers for maximum access speed
    const double* p_x = x.memptr();
    const arma::uword* p_ind = ind.memptr();
    double* p_out = out.memptr();
    arma::uword start = 0;
    for (arma::uword i = 0; i < m; ++i) {
        arma::uword end = p_ind[i];
        // Safety bound check to prevent crashing R if an index is too large
        if (end >= x.n_elem) end = x.n_elem - 1;
        double current_sum = 0.0;
        // Sum the elements for the current segment
        for (arma::uword j = start; j <= end; ++j) {
            current_sum += p_x[j];
        }
        p_out[i] = current_sum;
        // Next group starts exactly one element after the current group ends
        start = end + 1;
    }
    return out;
}

vec group_sum_old(const vec &x, const uvec &ind) {
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
      uword first_col = ind_RE.at(i).front();
      uword last_col  = ind_RE.at(i).back();
      out.at(i) = b.cols(first_col, last_col);
  }
  return out;
}

inline void mat2field_inplace (field<mat> &out, const mat &b, const field<uvec> &ind_RE) {
    uword n = ind_RE.n_elem;
    for (uword i = 0; i < n; ++i) {
        uword first_col = ind_RE.at(i).front();
        uword last_col  = ind_RE.at(i).back();
        out.at(i) = b.cols(first_col, last_col);
    }
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
                double tau, bool shrink) {
  vec z = x - mean;
  if (shrink) {
    Tau.diag() = lambda;
  }
  return -0.5 * tau * arma::as_scalar(z.t() * Tau * z);
}

vec propose_norm (const vec &thetas, const vec &scale, const uword &i) {
  vec proposed_thetas = thetas;
  proposed_thetas.at(i) = R::rnorm(thetas.at(i), scale.at(i));
  return proposed_thetas;
}

vec propose_unif (const vec &thetas, const vec &scale, uword i) {
  vec proposed_thetas = thetas;
  proposed_thetas.at(i) = R::runif(thetas.at(i) - Const_Unif_Proposal * scale.at(i),
                     thetas.at(i) + Const_Unif_Proposal * scale.at(i));
  return proposed_thetas;
}

vec propose_lnorm (const vec &thetas, const double &log_mu_i, const vec &scale,
                   uword i) {
  vec proposed_thetas = thetas;
  proposed_thetas.at(i) = R::rlnorm(log_mu_i, scale.at(i));
  return proposed_thetas;
}

vec propose_norm_mala (const vec &thetas, const vec &scale,
                       const double &deriv, uword i) {
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

vec propose_rnorm_vec (const vec &thetas, const vec &scale) {
  vec proposed_thetas = thetas;
  proposed_thetas = scale % randn(thetas.n_rows) + thetas;
  return proposed_thetas;
}

vec propose_mvnorm_vec (const mat &Sigma) {
  mat U = chol(Sigma, "lower");
  vec res = U * randn(U.n_cols);
  return res;
}

void mu_fun (arma::vec &eta, const std::string &link) {
    if (link == "identity") {
        return;
    } else if (link == "inverse") {
        eta = 1.0 / eta;
    } else if (link == "logit") {
        eta = 1.0 / (1.0 + arma::trunc_exp(-eta));
    } else if (link == "probit") {
        eta = arma::normcdf(eta);
    } else if (link == "cloglog") {
        eta = 1.0 - arma::trunc_exp(-arma::trunc_exp(eta));
    } else if (link == "log") {
        eta = arma::trunc_exp(eta);
    }
}

arma::vec log_dbinom (const arma::vec &x, const arma::vec &size, const arma::vec &prob) {
    arma::uword n = x.n_elem;
    arma::vec out(n, arma::fill::none);
    const double* px = x.memptr();
    const double* ps = size.memptr();
    const double* pp = prob.memptr();
    double* pout = out.memptr();
    for (arma::uword i = 0; i < n; ++i) {
        double xi = px[i];
        double ni = ps[i];
        double pi = pp[i];
        if (xi == 0.0) {
            pout[i] = ni * std::log(1.0 - pi);
        } else if (xi == ni) {
            pout[i] = ni * std::log(pi);
        } else {
            pout[i] = std::lgamma(ni + 1.0) - std::lgamma(xi + 1.0) - std::lgamma(ni - xi + 1.0)
            + xi * std::log(pi) + (ni - xi) * std::log(1.0 - pi);
        }
    }
    return out;
}

arma::vec log_dpois (const arma::vec &x, const arma::vec &lambda) {
    arma::uword n = x.n_elem;
    arma::vec out(n, arma::fill::none);
    const double* px = x.memptr();
    const double* pl = lambda.memptr();
    double* pout = out.memptr();
    for (arma::uword i = 0; i < n; ++i) {
        double xi = px[i];
        double lambda_i = pl[i];
        // Boundary short-circuit: log(Poisson(0; lambda)) = -lambda
        // Completely skips std::log() and std::lgamma() for zero counts!
        if (xi == 0.0) {
            pout[i] = -lambda_i;
        } else {
            // General case: x * log(lambda) - lambda - log(x!)
            pout[i] = xi * std::log(lambda_i) - lambda_i - std::lgamma(xi + 1.0);
        }
    }
    return out;
}

arma::vec log_dbbinom (const arma::vec &x, const arma::vec &size,
                       const arma::vec &prob, const double phi) {
    arma::uword n = x.n_elem;
    arma::vec out(n, arma::fill::none);
    const double* px = x.memptr();
    const double* ps = size.memptr();
    const double* pp = prob.memptr();
    double* pout = out.memptr();
    double lgamma_phi = std::lgamma(phi);
    for (arma::uword i = 0; i < n; ++i) {
        double xi = px[i];
        double ni = ps[i];
        double pi = pp[i];
        double A = phi * pi;
        double B = phi * (1.0 - pi);
        // Pre-compute the denominator for the beta function's numerator.
        // It only depends on 'ni' and 'phi', and is used in every branch.
        double lgamma_n_phi = std::lgamma(ni + phi);
        if (xi == 0.0) {
            // Short-circuit for x = 0
            // log_binom is 0. log(Gamma(A)) completely cancels out.
            pout[i] = std::lgamma(ni + B) - lgamma_n_phi + lgamma_phi - std::lgamma(B);
        } else if (xi == ni) {
            // Short-circuit for x = n
            // log_binom is 0. log(Gamma(B)) completely cancels out.
            pout[i] = std::lgamma(ni + A) - lgamma_n_phi + lgamma_phi - std::lgamma(A);
        } else {
            // General case for 0 < x < size
            // 1. Log Binomial Coefficient: log( n! / (x! * (n-x)!) )
            double log_binom = std::lgamma(ni + 1.0) - std::lgamma(xi + 1.0) - std::lgamma(ni - xi + 1.0);
            // 2. Log Beta Numerator: log( B(x+A, n-x+B) )
            double log_beta_num = std::lgamma(xi + A) + std::lgamma(ni - xi + B) - lgamma_n_phi;
            // 3. Log Beta Denominator: log( B(A, B) )
            double log_beta_den = std::lgamma(A) + std::lgamma(B) - lgamma_phi;
            pout[i] = log_binom + log_beta_num - log_beta_den;
        }
    }
    return out;
}

arma::vec log_dbernoulli (const arma::vec &x, const arma::vec &prob) {
    arma::uword n = x.n_elem;
    arma::vec out(n, arma::fill::none);
    const double* px = x.memptr();
    const double* pp = prob.memptr();
    double* pout = out.memptr();
    for (arma::uword i = 0; i < n; ++i) {
        double p_val = (px[i] > 0.0) ? pp[i] : (1.0 - pp[i]);
        pout[i] = std::log(p_val);
    }
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
    double constant = -std::log(sigma) - half_log2pi;
    double var2 = 2.0 * sigma * sigma;
    return constant - arma::square(x - mu) / var2;
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
    uword n = x.n_elem;
    vec out(n, fill::none);
    const double* px = x.memptr();
    const double* pscale = scale.memptr();
    double* pout = out.memptr();
    double lgamma_shape = std::lgamma(shape);
    double shape_m1 = shape - 1.0;
    for (uword i = 0; i < n; ++i) {
        double x_val = px[i];
        double scale_val = pscale[i];
        pout[i] = -lgamma_shape - shape * std::log(scale_val)
            + shape_m1 * std::log(x_val) - x_val / scale_val;
    }
    return out;
}

arma::vec log_dbeta (const arma::vec &x, const arma::vec &shape1, const arma::vec &shape2) {
    arma::uword n = x.n_elem;
    arma::vec out(n, arma::fill::none);
    const double* px = x.memptr();
    const double* psh1 = shape1.memptr();
    const double* psh2 = shape2.memptr();
    double* pout = out.memptr();
    for (arma::uword i = 0; i < n; ++i) {
        double xi = px[i];
        double a = psh1[i];
        double b = psh2[i];
        pout[i] = (a - 1.0) * std::log(xi) + (b - 1.0) * std::log(1.0 - xi) -
            std::lgamma(a) - std::lgamma(b) + std::lgamma(a + b);
    }
    return out;
}

vec log_dmvnrm_chol (const mat &x, const mat &L) {
    uword const n = x.n_rows, k = x.n_cols;
    vec out(n, arma::fill::none);
    mat V = inv(trimatu(L));
    double log_det = -arma::sum(arma::log(L.diag()));
    double constants = -(double)k / 2.0 * log2pi;
    double other_terms = constants + log_det;
    std::vector<double> z(k);
    for (uword i = 0; i < n; ++i) {
        for(uword j = 0; j < k; ++j) {
            z[j] = x.at(i, j);
        }
        for (uword j = k; j-- > 0;) {
            double tmp = 0.0;
            const double* col_j = V.colptr(j);
            for (uword c = 0; c <= j; ++c) {
                tmp += col_j[c] * z[c];
            }
            z[j] = tmp;
        }
        double sq_dist = 0.0;
        for(uword j = 0; j < k; ++j) {
            sq_dist += z[j] * z[j];
        }
        out[i] = other_terms - 0.5 * sq_dist;
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

inline void linpred_surv_inplace(field<mat> &eta,
                                 const field<mat> &X, const field<vec> &betas,
                                 const field<mat> &Z, const field<mat> &b,
                                 const uvec &id) {
    uword n_outcomes = X.n_elem;
    for (uword i = 0; i < n_outcomes; ++i) {
        const mat& X_i = X.at(i);
        const vec& betas_i = betas.at(i);
        const mat& Z_i = Z.at(i);
        const mat& b_i = b.at(i);
        uword n_betas = betas_i.n_rows;
        uword n_REs = b_i.n_cols;
        uword n_forms = X_i.n_cols / n_betas;
        uword N = X_i.n_rows;
        for (uword j = 0; j < n_forms; ++j) {
            uword j_beta_start = j * n_betas;
            uword j_beta_end = j_beta_start + n_betas - 1;
            eta.at(i).col(j) = X_i.cols(j_beta_start, j_beta_end) * betas_i;
            uword j_RE_start = j * n_REs;
            double* eta_col_ptr = eta.at(i).colptr(j);
            const uword* id_ptr = id.memptr();
            for (uword k = 0; k < n_REs; ++k) {
                const double* Z_col_ptr = Z_i.colptr(j_RE_start + k);
                const double* b_col_ptr = b_i.colptr(k);
                for (uword obs = 0; obs < N; ++obs) {
                    eta_col_ptr[obs] += Z_col_ptr[obs] * b_col_ptr[id_ptr[obs]];
                }
            }
        }
    }
}

field<mat> linpred_surv (const field<mat> &X, const field<vec> &betas,
                         const field<mat> &Z, const field<mat> &b,
                         const uvec &id) {
    uword n_outcomes = X.n_elem;
    field<mat> out(n_outcomes);
    for (uword i = 0; i < n_outcomes; ++i) {
        const mat& X_i = X.at(i);
        const vec& betas_i = betas.at(i);
        const mat& Z_i = Z.at(i);
        const mat& b_i = b.at(i);
        uword n_betas = betas_i.n_rows;
        uword n_REs = b_i.n_cols;
        uword n_forms = X_i.n_cols / n_betas;
        out.at(i).set_size(X_i.n_rows, n_forms);
        mat b_id = b_i.rows(id);
        for (uword j = 0; j < n_forms; ++j) {
            uword j_beta_start = j * n_betas;
            uword j_beta_end = j_beta_start + n_betas - 1;
            uword j_RE_start = j * n_REs;
            uword j_RE_end = j_RE_start + n_REs - 1;
            out.at(i).col(j) = X_i.cols(j_beta_start, j_beta_end) * betas_i +
                arma::sum(Z_i.cols(j_RE_start, j_RE_end) % b_id, 1);
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
    const mat& X_i = X.at(i);
    const vec& betas_i = betas.at(i);
    const mat& Z_i = Z.at(i);
    const mat& b_i = b.at(i);
    const uvec& id_i = id.at(i);
    out.at(i) = X_i * betas_i + arma::sum(Z_i % b_i.rows(id_i), 1);
  }
  return out;
}

inline void linpred_mixed_inplace (field<vec> &eta,
                                   const field<mat> &X,
                                   const field<vec> &betas,
                                   const field<mat> &Z,
                                   const field<mat> &b,
                                   const field<uvec> &id) {
    uword n_outcomes = X.n_elem;
    for (uword i = 0; i < n_outcomes; ++i) {
        const mat& X_i = X.at(i);
        const vec& betas_i = betas.at(i);
        const mat& Z_i = Z.at(i);
        const mat& b_i = b.at(i);
        const uvec& id_i = id.at(i);
        eta.at(i) = X_i * betas_i;
        uword N = Z_i.n_rows;
        uword q = Z_i.n_cols;
        double* eta_ptr = eta.at(i).memptr();
        const uword* id_ptr = id_i.memptr();
        for (uword k = 0; k < q; ++k) {
            const double* Z_col = Z_i.colptr(k);
            const double* b_col = b_i.colptr(k);
            for (uword obs = 0; obs < N; ++obs) {
                eta_ptr[obs] += Z_col[obs] * b_col[id_ptr[obs]];
            }
        }
    }
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

mat transf_eta (const mat &eta, const CharacterVector &fun_nams) {
    uword k = fun_nams.length();
    mat out(eta.n_rows, k, fill::none);
    for (uword i = 0; i < k; i++) {
        std::string fun = as<std::string>(fun_nams[i]);
        if (fun == "identity") {
            out.col(i) = eta;
        } else if (fun == "abs") {
            out.col(i) = abs(eta);
        } else if (fun == "expit") {
            out.col(i) = 1.0 / (1.0 + trunc_exp(-eta));
        } else if (fun == "exp" || fun == "dexp") {
            out.col(i) = trunc_exp(eta);
        } else if (fun == "dexpit") {
            mat pp = 1.0 / (1.0 + trunc_exp(-eta));
            out.col(i) = pp % (1.0 - pp);
        } else if (fun == "log") {
            out.col(i) = trunc_log(eta);
        } else if (fun == "log2") {
            out.col(i) = log2(eta);
        } else if (fun == "log10") {
            out.col(i) = log10(eta);
        } else if (fun == "sqrt") {
            out.col(i) = sqrt(eta);
        } else if (fun == "poly2") {
            out.col(i) = square(eta);
        } else if (fun == "poly3") {
            out.col(i) = eta % square(eta);
        } else if (fun == "poly4") {
            out.col(i) = square(square(eta));
        } else if (fun == "poly2(expit)") {
            out.col(i) = square(1.0 / (1.0 + trunc_exp(-eta)));
        } else if (fun == "poly3(expit)") {
            mat pp = 1.0 / (1.0 + trunc_exp(-eta));
            out.col(i) = pp % square(pp);
        } else if (fun == "poly4(expit)") {
            out.col(i) = square(square(1.0 / (1.0 + trunc_exp(-eta))));
        }
    }
    return out;
}

field<mat> create_Wlong(const field<mat> &eta, const field<mat> &U,
                        const field<uvec> &FunForms,
                        const List &Funs_FunForms) {
    uword n_outcomes = eta.n_elem;
    field<mat> out(n_outcomes);
    for (uword i = 0; i < n_outcomes; ++i) {
        out.at(i) = U.at(i);
        const mat& eta_i = eta.at(i);
        const uvec& FF_i = FunForms.at(i);
        List Funs_i = Funs_FunForms[i];
        uword current_col = 0;
        uword n_funs = Funs_i.length();
        for (uword j = 0; j < n_funs; ++j) {
            const vec& eta_col = eta_i.col(j);
            CharacterVector fun_nams = Funs_i[j];
            uword k = fun_nams.length();
            for (uword f = 0; f < k; ++f) {
                uword target_col = FF_i.at(current_col);
                std::string fun = as<std::string>(fun_nams[f]);
                auto out_col = out.at(i).col(target_col);
                if (fun == "identity") {
                    out_col %= eta_col;
                } else if (fun == "abs") {
                    out_col %= arma::abs(eta_col);
                } else if (fun == "expit") {
                    out_col %= 1.0 / (1.0 + trunc_exp(-eta_col));
                } else if (fun == "exp" || fun == "dexp") {
                    out_col %= trunc_exp(eta_col);
                } else if (fun == "dexpit") {
                    vec pp = 1.0 / (1.0 + trunc_exp(-eta_col));
                    out_col %= pp % (1.0 - pp);
                } else if (fun == "log") {
                    out_col %= trunc_log(eta_col);
                } else if (fun == "log2") {
                    out_col %= arma::log2(eta_col);
                } else if (fun == "log10") {
                    out_col %= arma::log10(eta_col);
                } else if (fun == "sqrt") {
                    out_col %= arma::sqrt(eta_col);
                } else if (fun == "poly2") {
                    out_col %= arma::square(eta_col);
                } else if (fun == "poly3") {
                    out_col %= eta_col % arma::square(eta_col);
                } else if (fun == "poly4") {
                    out_col %= arma::square(arma::square(eta_col));
                } else if (fun == "poly2(expit)") {
                    out_col %= arma::square(1.0 / (1.0 + trunc_exp(-eta_col)));
                } else if (fun == "poly3(expit)") {
                    vec pp = 1.0 / (1.0 + trunc_exp(-eta_col));
                    out_col %= pp % arma::square(pp);
                } else if (fun == "poly4(expit)") {
                    out_col %= arma::square(arma::square(1.0 / (1.0 + trunc_exp(-eta_col))));
                }
                current_col++;
            }
        }
    }
    return out;
}

mat calculate_Wlong_old (const field<mat> &X, const field<mat> &Z,
                     const field<mat> &U, const mat &Wlong_bar,
                     const mat &Wlong_sds,
                     const field<vec> &betas, const field<mat> &b,
                     const uvec &id, const field<uvec> &FunForms,
                     const List &Funs_FunForms) {
  field<mat> eta = linpred_surv(X, betas, Z, b, id);
  mat Wlong =
    docall_cbindF(create_Wlong(eta, U, FunForms, Funs_FunForms));
  Wlong.each_row() -= Wlong_bar;
  Wlong.each_row() /= Wlong_sds;
  return Wlong;
}

mat calculate_Wlong (const field<mat> &X, const field<mat> &Z,
                     const field<mat> &U, const mat &Wlong_bar,
                     const mat &Wlong_sds,
                     const field<vec> &betas, const field<mat> &b,
                     const uvec &id, const field<uvec> &FunForms,
                     const List &Funs_FunForms) {
    //field<mat> eta = linpred_surv(X, betas, Z, b, id);
    uword n_outcomes_X = X.n_elem;
    field<mat> eta(n_outcomes_X);
    for (uword i = 0; i < n_outcomes_X; ++i) {
        uword n_forms = X.at(i).n_cols / betas.at(i).n_rows;
        eta.at(i).set_size(X.at(i).n_rows, n_forms);
    }
    linpred_surv_inplace(eta, X, betas, Z, b, id);
    uword n_outcomes = U.n_elem;
    uword N = U.at(0).n_rows;
    uword total_cols = 0;
    uvec col_starts(n_outcomes);
    for (uword i = 0; i < n_outcomes; ++i) {
        col_starts.at(i) = total_cols;
        total_cols += U.at(i).n_cols;
    }
    mat Wlong(N, total_cols, arma::fill::none);
    for (uword i = 0; i < n_outcomes; ++i) {
        uword start_col = col_starts.at(i);
        uword end_col = start_col + U.at(i).n_cols - 1;
        Wlong.cols(start_col, end_col) = U.at(i);
        const mat& eta_i = eta.at(i);
        const uvec& FF_i = FunForms.at(i);
        List Funs_i = Funs_FunForms[i];
        uword current_col = 0;
        uword n_funs = Funs_i.length();
        for (uword j = 0; j < n_funs; ++j) {
            const vec& eta_col = eta_i.col(j);
            IntegerVector fun_ints = Funs_i[j];
            uword k = fun_ints.length();
            for (uword f = 0; f < k; ++f) {
                uword target_col = start_col + FF_i.at(current_col);
                auto out_col = Wlong.col(target_col);
                int fun = fun_ints[f];
                switch (fun) {
                case 1: // "identity"
                    out_col %= eta_col;
                    break;
                case 2: // "abs"
                    out_col %= arma::abs(eta_col);
                    break;
                case 3: // "expit"
                    out_col %= 1.0 / (1.0 + trunc_exp(-eta_col));
                    break;
                case 4: // "exp" or "dexp"
                    out_col %= trunc_exp(eta_col);
                    break;
                case 5: // "dexpit"
                    { // Brackets keep 'pp' safely scoped within this case
                        vec pp = 1.0 / (1.0 + trunc_exp(-eta_col));
                        out_col %= pp % (1.0 - pp);
                        break;
                    }
                case 6: // "log"
                    out_col %= trunc_log(eta_col);
                    break;
                case 7: // "log2"
                    out_col %= arma::log2(eta_col);
                    break;
                case 8: // "log10"
                    out_col %= arma::log10(eta_col);
                    break;
                case 9: // "sqrt"
                    out_col %= arma::sqrt(eta_col);
                    break;
                case 10: // "poly2"
                    out_col %= arma::square(eta_col);
                    break;
                case 11: // "poly3"
                    out_col %= eta_col % arma::square(eta_col);
                    break;
                case 12: // "poly4"
                    out_col %= arma::square(arma::square(eta_col));
                    break;
                case 13: // "poly2(expit)"
                    out_col %= arma::square(1.0 / (1.0 + trunc_exp(-eta_col)));
                    break;
                case 14: // "poly3(expit)"
                    {
                        vec pp = 1.0 / (1.0 + trunc_exp(-eta_col));
                        out_col %= pp % arma::square(pp);
                        break;
                    }
                case 15: // "poly4(expit)"
                    out_col %= arma::square(arma::square(1.0 / (1.0 + trunc_exp(-eta_col))));
                    break;
                default:
                    Rcpp::stop("Unknown transformation function integer code.");
                }
                current_col++;
            }
        }
    }
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

mat rank1_update_old (const mat &U, // performs rank-1 update. If U = chol(M), returns chol(M + v * v.t())
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

mat rank1_update (const mat &U, const vec &v) {
    uword n = v.n_elem;
    mat Res = U;
    vec v2 = v;
    for (uword i = 0; i < n; ++i) {
        double res_ii = Res.at(i, i);
        double v2_i = v2.at(i);
        double r = std::sqrt(res_ii * res_ii + v2_i * v2_i);
        double c = r / res_ii;
        double s = v2_i / res_ii;
        Res.at(i, i) = r;
        for (uword j = i + 1; j < n; ++j) {
            double res_ij = Res.at(i, j);
            double v2_j   = v2.at(j);
            double new_res_ij = (res_ij + s * v2_j) / c;
            Res.at(i, j) = new_res_ij;
            v2.at(j) = c * v2_j - s * new_res_ij;
        }
    }
    return Res;
}

mat chol_update_old(const mat &U, // If U = chol(M), returns chol(M.submat(keep, keep))
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

mat chol_update (const mat &U, const uvec &keep) {
    uword N = U.n_cols;
    mat Res = U; // Only one deep copy made
    // 1. Generate the removal indices exactly once
    uvec rem = regspace<uvec>(0, N - 1);
    rem.shed_rows(keep);
    // 2. Process each removal index without ever shrinking the matrix
    for (uword i = 0; i < rem.n_elem; ++i) {
        uword k = rem.at(i);
        // 3. INLINED RANK-1 UPDATE
        // We map the submatrix directly to the global coordinates of Res.
        for (uword row_idx = k + 1; row_idx < N; ++row_idx) {
            double res_ii = Res.at(row_idx, row_idx);
            // Res.at(k, row_idx) perfectly acts as our 'v' vector!
            double v_i    = Res.at(k, row_idx);
            // Fast hypotenuse math (avoids pow)
            double r = std::sqrt(res_ii * res_ii + v_i * v_i);
            double c = r / res_ii;
            double s = v_i / res_ii;
            Res.at(row_idx, row_idx) = r;
            for (uword col_idx = row_idx + 1; col_idx < N; ++col_idx) {
                double res_ij = Res.at(row_idx, col_idx);
                double v_j    = Res.at(k, col_idx);
                double new_res_ij = (res_ij + s * v_j) / c;
                Res.at(row_idx, col_idx) = new_res_ij;
                // Overwrite the removed row to act as our updated 'v' for the next step
                Res.at(k, col_idx) = c * v_j - s * new_res_ij;
            }
        }
    }
    // 4. Extract the final submatrix in one single allocation
    return Res.submat(keep, keep);
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

arma::uvec std_setdiff(arma::uvec& x, arma::uvec& y) {
    std::vector<int> a = arma::conv_to< std::vector<int> >::from(arma::sort(x));
    std::vector<int> b = arma::conv_to< std::vector<int> >::from(arma::sort(y));
    std::vector<int> out;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(out, out.end()));
    return arma::conv_to<arma::uvec>::from(out);
}

double logSumExp (const vec &x) {
    double maxval = max(x);
    double out = maxval + log(sum(exp(x - maxval)));
    return out;
}

vec lse (const vec &xx, const uvec &id_h2, const uvec & intgr_ind) {
    uvec unq_idh2 = unique(id_h2);
    uword n = unq_idh2.n_rows;
    vec out(n);
    for (uword i = 0; i < n; ++i) {
        uvec idx = find(id_h2 == i);
        uvec intgr_i = intgr_ind(idx);
        uword nn = intgr_i.n_rows;
        if (nn > 0) {
            out(i) = logSumExp(xx(idx));
        } else {
            out(i) = sum(xx(idx));
        }
    }
    return out;
}

#endif

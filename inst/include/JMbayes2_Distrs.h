#ifndef JMBAYES2DISTRS
#define JMBAYES2DISTRS

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

static double const log2pi = std::log(2.0 * M_PI);
static double const Const_Unif_Proposal = 0.5 * std::pow(12.0, 0.5);

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
                          const int& k, const int& i) {
  field<vec> proposed_thetas = thetas;
  proposed_thetas.at(k).at(i) = R::rnorm(thetas.at(k).at(i),
                     scale.at(k).at(i));
  return proposed_thetas;
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

#endif

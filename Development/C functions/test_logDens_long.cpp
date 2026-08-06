#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec lbeta_arma_fast (const vec &a, const vec &b) {
    uword n = a.n_elem;
    vec out(n, fill::none);
    const double* pa = a.memptr();
    const double* pb = b.memptr();
    double* pout = out.memptr();
    for (uword i = 0; i < n; ++i) {
        pout[i] = R::lbeta(pa[i], pb[i]);
    }
    return out;
}

// [[Rcpp::export]]
vec lchoose_arma_fast (const vec &n, const vec &k) {
    uword n_elem = n.n_elem;
    vec out(n_elem, fill::none);
    const double* pn = n.memptr();
    const double* pk = k.memptr();
    double* pout = out.memptr();
    for (uword i = 0; i < n_elem; ++i) {
        pout[i] = R::lchoose(pn[i], pk[i]);
    }
    return out;
}

// [[Rcpp::export]]
vec log_dbinom_fast (const vec &x, const vec &size, const vec &prob) {
    uword n = x.n_elem;
    vec out(n, fill::none);
    const double* px = x.memptr();
    const double* ps = size.memptr();
    const double* pp = prob.memptr();
    double* pout = out.memptr();
    for (uword i = 0; i < n; ++i) {
        pout[i] = R::dbinom(px[i], ps[i], pp[i], 1);
    }
    return out;
}

// [[Rcpp::export]]
arma::vec log_dbinom_ultra (const arma::vec &x, const arma::vec &size, const arma::vec &prob) {
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

        // Boundary short-circuits: handles x = 0 and x = size without computing gamma functions
        if (xi == 0.0) {
            pout[i] = ni * std::log(1.0 - pi);
        } else if (xi == ni) {
            pout[i] = ni * std::log(pi);
        } else {
            // General case for 0 < x < size
            pout[i] = std::lgamma(ni + 1.0) - std::lgamma(xi + 1.0) - std::lgamma(ni - xi + 1.0)
            + xi * std::log(pi) + (ni - xi) * std::log(1.0 - pi);
        }
    }

    return out;
}

// [[Rcpp::export]]
vec log_dpois_fast (const vec &x, const vec &lambda) {
    uword n = x.n_elem;
    vec out(n, fill::none);
    const double* px = x.memptr();
    const double* pl = lambda.memptr();
    double* pout = out.memptr();
    for (uword i = 0; i < n; ++i) {
        pout[i] = R::dpois(px[i], pl[i], 1);
    }
    return out;
}

// 2. Loop Fusion: Eliminate intermediate vector allocations
// [[Rcpp::export]]
vec log_dbbinom_fast (const vec &x, const vec &size, const vec &prob, const double &phi) {
    uword n = x.n_elem;
    vec out(n, fill::none);
    const double* px = x.memptr();
    const double* ps = size.memptr();
    const double* pp = prob.memptr();
    double* pout = out.memptr();

    // Compute element-wise in a single pass without allocating A, B, log_numerator, etc.
    for (uword i = 0; i < n; ++i) {
        double A = phi * pp[i];
        double B = phi * (1.0 - pp[i]);
        pout[i] = R::lchoose(ps[i], px[i]) +
            R::lbeta(px[i] + A, ps[i] - px[i] + B) -
            R::lbeta(A, B);
    }
    return out;
}

// [[Rcpp::export]]
arma::vec log_dbbinom_ultra(const arma::vec &x,
                            const arma::vec &size,
                            const arma::vec &prob,
                            const double phi) {
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

// [[Rcpp::export]]
vec log_dnbinom_fast (const vec &x, const vec &mu, const double &size) {
    uword n = x.n_elem;
    vec out(n, fill::none);
    const double* px = x.memptr();
    const double* pmu = mu.memptr();
    double* pout = out.memptr();

    // Precompute constant scalar values outside the loop!
    double lgamma_size = std::lgamma(size);
    double size_log_size = size * std::log(size);

    // Single pass calculation
    for (uword i = 0; i < n; ++i) {
        double mu_val = pmu[i];
        double x_val = px[i];

        double log_mu_size = std::log(mu_val + size);
        double comp1 = std::lgamma(x_val + size) - lgamma_size - std::lgamma(x_val + 1.0);
        double comp2 = size_log_size - size * log_mu_size;
        double comp3 = x_val * (std::log(mu_val) - log_mu_size);

        pout[i] = comp1 + comp2 + comp3;
    }
    return out;
}

// 3. Mathematical simplification and Expression Templates
// [[Rcpp::export]]
vec log_dnorm_fast (const vec &x, const vec &mu, const double &sigma) {
    // Avoid creating a vector of 'sigmas'.
    // The log-PDF of a normal distribution is mathematically:
    // -0.5 * log(2 * pi) - log(sigma) - (x - mu)^2 / (2 * sigma^2)

    double constant = -std::log(sigma) - 0.5 * std::log(2.0 * M_PI);
    double var2 = 2.0 * sigma * sigma;

    // Armadillo computes this directly into the output vector with 0 temporaries
    return constant - arma::square(x - mu) / var2;
}

// [[Rcpp::export]]
vec lbeta_arma (const vec &a, const vec &b) {
    uword n = a.n_rows;
    vec out(n);
    for (uword i = 0; i < n; ++i) {
        out.at(i) = R::lbeta(a.at(i), b.at(i));
    }
    return out;
}

// [[Rcpp::export]]
vec lchoose_arma (const vec &n, const vec &k) {
    uword n_ = n.n_rows;
    vec out(n_);
    for (uword i = 0; i < n_; ++i) {
        out.at(i) = R::lchoose(n.at(i), k.at(i));
    }
    return out;
}

// [[Rcpp::export]]
vec log_dbinom (const vec &x, const vec &size, const vec &prob) {
    uword n = x.n_rows;
    vec out(n);
    for (uword i = 0; i < n; ++i) {
        out.at(i) = R::dbinom(x.at(i), size.at(i), prob.at(i), 1);
    }
    return out;
}

// [[Rcpp::export]]
vec log_dpois (const vec &x, const vec &lambda) {
    uword n = x.n_rows;
    vec out(n);
    for (uword i = 0; i < n; ++i) {
        out.at(i) = R::dpois(x.at(i), lambda.at(i), 1);
    }
    return out;
}

// [[Rcpp::export]]
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

// [[Rcpp::export]]
vec log_dnbinom (const vec &x, const vec &mu, const double &size) {
    vec log_mu_size = log(mu + size);
    vec comp1 = lgamma(x + size) - lgamma(size) - lgamma(x + 1);
    vec comp2 = size * log(size) - size * log_mu_size;
    vec comp3 = x % (log(mu) - log_mu_size);
    vec out = comp1 + comp2 + comp3;
    return out;
}

// [[Rcpp::export]]
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

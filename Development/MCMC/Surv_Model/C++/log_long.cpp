#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

vec mu_fun (const vec &eta, const std::string &link) {
    uword n = eta.n_rows;
    vec exp_eta(n);
    vec out(n);
    if (link == "identity") {
        out = eta;
    } else if (link == "logit") {
        exp_eta = trunc_exp(eta);
        out = exp_eta / (1 + exp_eta);
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

// [[Rcpp::export]]
vec log_long (const field<mat> &y, const field<vec> &eta, const vec &scales,
              const vec &extra_parms, const CharacterVector &families,
              const CharacterVector &links, const field<uvec> &ids) {
    uword n_outcomes = y.size();
    uword n = y.at(0).n_rows;
    vec log_contr(n);
    vec out(n);
    for (uword i = 0; i < n_outcomes; ++i) {
        uvec id_i = ids.at(i);
        mat y_i = y.at(i);
        vec mu_i = mu_fun(eta.at(i), links[i]);
        double scale_i = scales.at(i);
        double extr_prm_i = extra_parms.at(i);
        if (families[i] == "gaussian") {
            log_contr = log_dnorm(y_i, mu_i, scale_i);
        } else if (families[i] == "binomial") {
            uword k = y_i.n_cols;
            if (k == 2) {
                log_contr = log_dbinom(y.col(0), y.col(0) + y.col(1), mu_i);
            } else {
                log_contr = y_i % log(mu_i) + (1 - y_i) % log(1 - mu_i);
            }
        } else if (families[i] == "poisson") {
            log_contr = log_dpois(y_i, mu_i);
        } else if (families[i] == "negative binomial") {
            log_contr = log_dnbinom(y_i, mu_i, scale_i);
        } else if (families[i] == "beta") {
            log_contr = log_dbeta(y_i, mu_i * scale_i, scale_i * (1 - mu_i));
        } else if (families[i] == "Student-t") {
            log_contr = log_dt((y_i - mu_i) / scale_i, extr_prm_i) - log(scale_i);
        } else if (families[i] == "Gamma") {
            log_contr = log_dgamma(y_i, square(mu_i) / scale_i, scale_i / mu_i);
        } else if (families[i] == "beta binomial") {
            uword k = y_i.n_cols;
            if (k == 2) {
                log_contr = log_dbbinom(y_i.col(0), y_i.col(0) + y_i.col(1), mu_i, scale_i);
            } else {
                vec ones(n, fill::ones);
                log_contr = log_dbbinom(y_i, ones, mu_i, scale_i);
            }
        } else if (families[i] == "unit Lindley") {
            vec theta = 1 / mu_i - 1;
            vec comp1 = 2 * log(theta) - log(1 + theta);
            vec comp2 = - 3 * log(1 - y_i);
            vec comp3 = - (theta * y_i) / (1 - y_i);
            log_contr = comp1 + comp2 + comp3;
        }
        out += group_sum(log_contr, id_i)
    }
    return out;
}

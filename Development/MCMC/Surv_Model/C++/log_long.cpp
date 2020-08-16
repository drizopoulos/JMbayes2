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

// [[Rcpp::export]]
vec log_dt_arma (const vec &x, const double &df) {
    uword n = x.n_rows;
    vec out(n);
    for (uword i = 0; i < n; ++i) {
        out.at(i) = R::dt(x.at(i), df, 1);
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

/*
vec log_long (const field<mat> &y, const field<vec> &eta, const vec &scales, const vec &extra_parms,
              const CharacterVector &families, const CharacterVector &links, const field<uvec> &ids) {
    uword n_outcomes = y.size();
    uword n = y.at(0).n_rows;
    vec log_contr(n);
    vec out(n);
    for (uword i = 0; i < n_outcomes; ++i) {
        uvec id_i = id.at(i);
        mat y_i = y.at(i);
        vec mu_i = mu_fun(eta.at(i), links[i]);
        double scale_i = scales.at(i);
        double extr_prm_i = extra_parms.at(i);
        if (families[i] == "gaussian") {
            log_contr = log_normpdf(y_i, mu_i, scale_i);
        } else if (families[i] == "binomial") {
            uword k = y_i.n_cols;
            if (k == ) {
                log_contr = dbinom(y.col(0), y.col(0) + y.col(1), mu_i, true); // do y.col(0) + y.col(1) once
            } else {
                log_contr = y_i % log(mu_i) + (1 - y_i) % log(1 - mu_i);
            }
        } else if (families[i] == "poisson") {
            log_contr = dpois(y_i, mu_i, true);
        } else if (families[i] == "negative binomial") {
            log_contr = dnbinom_mu(y_i, scale_i, mu_i, true);
        } else if (families[i] == "beta") {
            log_contr = dbeta(y_i, mu_i * scale_i, scale_i * (1 - mu_i), true);
        } else if (families[i] == "Student-t") {
            log_contr = dt((y_i - mu_i) / scale_i, extr_prm_i, true) - log(scale_i);
        } else if (families[i] == "Gamma") {
            log_contr = dgamma(y_i, square(mu_i) / scale_i, scale_i / mu_i, true);
        } else if (families[i] == "beta binomial") {
            uword k = y_i.n_cols;
            if (k == 2) {
                log_contr = log_dbbinom(y_i.col(0), y_i.col(0) + y_i.col(1), mu_i, scale_i);
            } else {
                log_contr = log_dbbinom(y_i, 1 + 0 * y_i, mu_i, scale_i);
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
*/

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec mu_fun_fast (const arma::vec &eta, const std::string &link) {
    if (link == "identity") {
        return eta;
    } else if (link == "inverse") {
        return 1.0 / eta;
    } else if (link == "logit") {
        return 1.0 / (1.0 + arma::trunc_exp(-eta));
    } else if (link == "probit") {
        return arma::normcdf(eta);
    } else if (link == "cloglog") {
        return 1.0 - arma::trunc_exp(-arma::trunc_exp(eta));
    } else if (link == "log") {
        return arma::trunc_exp(eta);
    }
}

// [[Rcpp::export]]
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




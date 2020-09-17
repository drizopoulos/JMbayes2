#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
vec mu_fun (const vec &eta, const CharacterVector &link) {
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






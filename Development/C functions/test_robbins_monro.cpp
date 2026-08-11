#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double robbins_monro_fast(double scale, double acceptance_it,
                          int it, double target_acceptance = 0.45) {
    double it_d = static_cast<double>(it);
    if (acceptance_it > 0.0) {
        return scale * (1.0 + 1.0 / (target_acceptance * it_d));
    } else {
        return scale * (1.0 - 1.0 / ((1.0 - target_acceptance) * it_d));
    }
}

// [[Rcpp::export]]
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

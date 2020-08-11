#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

static double const log2pi = std::log(2.0 * M_PI);
static double const Const_Unif_Proposal = 0.5 * std::pow(12.0, 0.5);

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

double logPC_D_sds (const vec &sds, const mat &L, const mat &b,
                    const double &prior_D_sds_df,
                    const double &prior_D_sds_sigma) {
    mat chol_Sigma = L.each_row() % sds.t();
    double log_p_b = sum(log_dmvnrm_chol(b, chol_Sigma));
    double log_p_sds = sum(log_dht(sds, prior_D_sds_sigma, prior_D_sds_df));
    double out = log_p_b + log_p_sds;
    return out;
}

double logPC_D_L (const mat &L, const vec &sds, const mat &b,
                  const double &prior_D_L_etaLKJ) {
    uword p = L.n_rows;
    mat chol_Sigma = L.each_row() % sds.t(); // check this
    double log_p_b = sum(log_dmvnrm_chol(b, chol_Sigma));
    double log_p_L(0.0);
    for (unsigned i = 1; i < p; ++i) {
        log_p_L += (p - i - 1.0 + 2.0 * prior_D_L_etaLKJ - 2.0) * log(L.at(i, i));
    }
    double out = log_p_b + log_p_L;
    return out;
}

vec propose_norm (const vec &thetas, const vec &scale, const uword &i) {
    // change R::rnorm to Armadillo native, and set the seed in Armadillo.
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

void update_D (mat &L, vec &sds, const mat &b,
               const uvec &upper_part,
               const double &prior_D_sds_df,
               const double &prior_D_sds_sigma,
               const double &prior_D_L_etaLKJ,
               const int &it, const bool &MALA,
               mat &res_sds, mat &res_L,
               vec &scale_sds, vec &scale_L,
               mat &acceptance_sds, mat &acceptance_L) {
    uword n_sds = sds.n_rows;
    uword n_L = upper_part.n_rows;
    double denominator_sds = logPC_D_sds(sds, L, b, prior_D_sds_df,
                                         prior_D_sds_sigma);
    for (uword i = 0; i < n_sds; ++i) {
        double SS = 0.5 * pow(scale_sds.at(i), 2.0);
        double log_mu_current = log(sds.at(i)) - SS;
        vec proposed_sds = propose_lnorm(sds, log_mu_current, scale_sds, i);
        double numerator_sds = logPC_D_sds(proposed_sds, L, b,
                                           prior_D_sds_df, prior_D_sds_sigma);
        double log_mu_proposed = log(proposed_sds.at(i)) - SS;
        double log_ratio_sds = numerator_sds - denominator_sds +
            R::dlnorm(sds.at(i), log_mu_proposed, scale_sds.at(i), true) -
            R::dlnorm(proposed_sds.at(i), log_mu_current, scale_sds.at(i), true);
        if (std::isfinite(log_ratio_sds) && exp(log_ratio_sds) > R::runif(0.0, 1.0)) {
            sds = proposed_sds;
            denominator_sds = numerator_sds;
            acceptance_sds.at(it, i) = 1;
        }
        if (it > 19) {
            scale_sds.at(i) =
                robbins_monro(scale_sds.at(i), acceptance_sds.at(it, i),
                              it);
        }
        res_sds.at(it, i) = sds.at(i);
    }
    double denominator_L = logPC_D_L(L, sds, b, prior_D_L_etaLKJ);
    for (uword i = 0; i < n_L; ++i) {
        uword upper_part_i = upper_part.at(i);
        double deriv_current(0.0);
        double mu_current(0.0);
        mat proposed_L = L;
        if (MALA) {
            deriv_current = deriv_L(L, sds, b, denominator_L, i, upper_part,
                                    prior_D_L_etaLKJ);
            mu_current = L.at(upper_part_i) + 0.5 * scale_L.at(i) * deriv_current;
            proposed_L = propose_L(L, scale_L, upper_part, deriv_current, i, true);
        } else {
            proposed_L = propose_L(L, scale_L, upper_part, deriv_current, i);
        }
        double numerator_L(0.0);
        double deriv_proposed(0.0);
        double mu_proposed(0.0);
        double log_ratio_L(0.0);
        if (proposed_L.is_finite()) {
            numerator_L = logPC_D_L(proposed_L, sds, b, prior_D_L_etaLKJ);
            if (MALA) {
                deriv_proposed = deriv_L(proposed_L, sds, b, numerator_L,
                                         i, upper_part, prior_D_L_etaLKJ);
                mu_proposed = proposed_L.at(upper_part_i) +
                    0.5 * scale_L.at(i) * deriv_proposed;
                log_ratio_L = numerator_L - denominator_L +
                    log_normpdf(L.at(upper_part_i), mu_proposed, sqrt(scale_L.at(i))) -
                    log_normpdf(proposed_L.at(upper_part_i), mu_current, sqrt(scale_L.at(i)));
            } else {
                log_ratio_L = numerator_L - denominator_L;
            }
        }
        if (std::isfinite(log_ratio_L) && exp(log_ratio_L) > R::runif(0.0, 1.0)) {
            L = proposed_L;
            denominator_L = numerator_L;
            acceptance_L.at(it, i) = 1;
        }
        if (it > 19) {
            scale_L.at(i) =
                robbins_monro(scale_L.at(i), acceptance_L.at(it, i),
                              it);
        }
        res_L.at(it, i) = L.at(upper_part_i);
    }
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
    int n_list = Mats.size();
    field<mat> res(n_list);
    for (int i = 0; i < n_list; ++i) {
        res.at(i) = as<mat>(Mats[i]);
    }
    return res;
}

field<vec> List2Field_vec (const List &Vecs) {
    int n_list = Vecs.size();
    field<vec> res(n_list);
    for (int i = 0; i < n_list; ++i) {
        res.at(i) = as<vec>(Vecs[i]);
    }
    return res;
}

field<mat> create_storage(const field<vec> &F, const int &n_iter) {
    int n = F.size();
    field<mat> out(n);
    for (int i = 0; i < n; ++i) {
        vec aa = F.at(i);
        int n_i = aa.n_rows;
        mat tt(n_iter, n_i, fill::zeros);
        out.at(i) = tt;
    }
    return out;
}

vec Wlong_alphas_fun(const field<mat> &Mats, const field<vec> &coefs) {
    int n = Mats.size();
    int n_rows = Mats.at(0).n_rows;
    vec out(n_rows, fill::zeros);
    for (int k = 0; k < n; ++k) {
        out += Mats.at(k) * coefs.at(k);
    }
    return out;
}

mat docall_cbind (const List &Mats_) {
    field<mat> Mats = List2Field_mat(Mats_);
    int n = Mats.size();
    uvec ncols(n);
    for (int k = 0; k < n; ++k) {
        ncols.at(k) = Mats.at(k).n_cols;
    }
    int N = sum(ncols);
    int col_start = 0;
    int col_end = ncols.at(0) - 1;
    mat out(Mats.at(0).n_rows, N);
    for (int k = 0; k < n; ++k) {
        if (k > 0) {
            col_start += ncols.at(k - 1);
            col_end += ncols.at(k);
        }
        out.cols(col_start, col_end) = Mats.at(k);
    }
    return out;
}

uvec create_fast_ind (const uvec &group) {
    unsigned int l = group.n_rows;
    uvec ind = find(group.rows(1, l - 1) != group.rows(0, l - 2));
    unsigned int k = ind.n_rows;
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

double log_density_surv (const vec &W0H_bs_gammas,
                         const vec &W0h_bs_gammas,
                         const vec &W0H2_bs_gammas,
                         const vec &WH_gammas,
                         const vec &Wh_gammas,
                         const vec &WH2_gammas,
                         const vec &WlongH_alphas,
                         const vec &Wlongh_alphas,
                         const vec &WlongH2_alphas,
                         const vec &log_Pwk, const vec &log_Pwk2,
                         const uvec &indFast_H,
                         const uvec &which_event,
                         const uvec &which_right_event,
                         const uvec &which_left,
                         const bool &any_interval,
                         const uvec &which_interval) {
    vec lambda_H = W0H_bs_gammas + WH_gammas + WlongH_alphas;
    vec H = group_sum(exp(log_Pwk + lambda_H), indFast_H);
    int n = H.n_rows;
    vec lambda_h(n);
    lambda_h.elem(which_event) = W0h_bs_gammas.elem(which_event) +
        Wh_gammas.elem(which_event) + Wlongh_alphas.elem(which_event);
    vec log_Lik_surv(n);
    log_Lik_surv.elem(which_right_event) = - H.elem(which_right_event);
    log_Lik_surv.elem(which_event) += lambda_h.elem(which_event);
    log_Lik_surv.elem(which_left) = log1p(- exp(- H.elem(which_left)));
    vec lambda_H2(lambda_H.n_rows);
    vec H2(n);
    if (any_interval) {
        lambda_H2 = W0H2_bs_gammas + WH2_gammas + WlongH2_alphas;
        H2 = group_sum(exp(log_Pwk2 + lambda_H2), indFast_H);
        log_Lik_surv.elem(which_interval) = - H.elem(which_interval) +
            log(- expm1(- H2.elem(which_interval)));
    }
    double logLik = sum(log_Lik_surv);
    return logLik;
}

void update_bs_gammas (vec& bs_gammas, vec& gammas, vec& alphas,
                       vec& W0H_bs_gammas, vec& W0h_bs_gammas, vec& W0H2_bs_gammas,
                       vec& WH_gammas, vec& Wh_gammas, vec& WH2_gammas,
                       vec& WlongH_alphas, vec& Wlongh_alphas, vec& WlongH2_alphas,
                       vec& log_Pwk, vec& log_Pwk2, uvec& id_H,
                       uvec& which_event, uvec& which_right_event, uvec& which_left, uvec& which_interval,
                       bool& any_event, bool& any_interval,
                       vec& prior_mean_bs_gammas, mat& prior_Tau_bs_gammas, double& tau_bs_gammas,
                       vec& prior_mean_gammas, mat& prior_Tau_gammas,
                       vec& prior_mean_alphas, mat& prior_Tau_alphas,
                       double& denominator_surv, int& it,
                       /////
                       mat& W0_H, mat& W0_h, mat& W0_H2,
                       vec& scale_bs_gammas,
                       mat& acceptance_bs_gammas,
                       mat& res_bs_gammas) {
    for (unsigned int i = 0; i < bs_gammas.n_rows; ++i) {
        vec proposed_bs_gammas = propose_norm(bs_gammas, scale_bs_gammas, i);
        vec proposed_W0H_bs_gammas = W0_H * proposed_bs_gammas;
        vec proposed_W0h_bs_gammas(W0_h.n_rows);
        vec proposed_W0H2_bs_gammas(W0_H2.n_rows);
        if (any_event) {
            proposed_W0h_bs_gammas = W0_h * proposed_bs_gammas;
        }
        if (any_interval) {
            proposed_W0H2_bs_gammas = W0_H2 * proposed_bs_gammas;
        }
        double numerator_surv =
            log_density_surv(proposed_W0H_bs_gammas, proposed_W0h_bs_gammas, proposed_W0H2_bs_gammas,
                             WH_gammas, Wh_gammas, WH2_gammas,
                             WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                             log_Pwk, log_Pwk2, id_H,
                             which_event, which_right_event, which_left,
                             any_interval, which_interval) +
                    logPrior(proposed_bs_gammas, prior_mean_bs_gammas,
                                          prior_Tau_bs_gammas, tau_bs_gammas) +
                    logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
                    logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
        double log_ratio = numerator_surv - denominator_surv;
        if (is_finite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
            bs_gammas = proposed_bs_gammas;
            W0H_bs_gammas = proposed_W0H_bs_gammas;
            if (any_event) {
                W0h_bs_gammas = proposed_W0h_bs_gammas;
            }
            if (any_interval) {
                W0H2_bs_gammas = proposed_W0H2_bs_gammas;
            }
            denominator_surv = numerator_surv;
            acceptance_bs_gammas.at(it, i) = 1;
        }
        if (it > 19) {
            scale_bs_gammas.at(i) =
                robbins_monro(scale_bs_gammas.at(i),
                              acceptance_bs_gammas.at(it, i), it);
        }
        res_bs_gammas.at(it, i) = bs_gammas.at(i);
    }
}

void update_gammas (vec& bs_gammas, vec& gammas, vec& alphas,
                    vec& W0H_bs_gammas, vec& W0h_bs_gammas, vec& W0H2_bs_gammas,
                    vec& WH_gammas, vec& Wh_gammas, vec& WH2_gammas,
                    vec& WlongH_alphas, vec& Wlongh_alphas, vec& WlongH2_alphas,
                    vec& log_Pwk, vec& log_Pwk2, uvec& id_H,
                    uvec& which_event, uvec& which_right_event, uvec& which_left, uvec& which_interval,
                    bool& any_event, bool& any_interval,
                    vec& prior_mean_bs_gammas, mat& prior_Tau_bs_gammas, double& tau_bs_gammas,
                    vec& prior_mean_gammas, mat& prior_Tau_gammas,
                    vec& prior_mean_alphas, mat& prior_Tau_alphas,
                    double& denominator_surv, int& it,
                    /////
                    mat& W_H, mat& W_h, mat& W_H2,
                    vec& scale_gammas,
                    mat& acceptance_gammas,
                    mat& res_gammas) {
    for (unsigned int i = 0; i < gammas.n_rows; ++i) {
        vec proposed_gammas = propose_norm(gammas, scale_gammas, i);
        vec proposed_WH_gammas = W_H * proposed_gammas;
        vec proposed_Wh_gammas(W_h.n_rows);
        vec proposed_WH2_gammas(W_H2.n_rows);
        if (any_event) {
            proposed_Wh_gammas = W_h * proposed_gammas;
        }
        if (any_interval) {
            proposed_WH2_gammas = W_H2 * proposed_gammas;
        }
        double numerator_surv =
            log_density_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                             proposed_WH_gammas, proposed_Wh_gammas, proposed_WH2_gammas,
                             WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                             log_Pwk, log_Pwk2, id_H,
                             which_event, which_right_event, which_left,
                             any_interval, which_interval) +
            logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas,
                     tau_bs_gammas) +
            logPrior(proposed_gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
            logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
        double log_ratio = numerator_surv - denominator_surv;
        if (is_finite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
            gammas = proposed_gammas;
            WH_gammas = proposed_WH_gammas;
            if (any_event) {
                Wh_gammas = proposed_Wh_gammas;
            }
            if (any_interval) {
                WH2_gammas = proposed_WH2_gammas;
            }
            denominator_surv = numerator_surv;
            acceptance_gammas.at(it, i) = 1;
        }
        if (it > 19) {
            scale_gammas.at(i) =
                robbins_monro(scale_gammas.at(i),
                              acceptance_gammas.at(it, i), it);
        }
        // store results
        res_gammas.at(it, i) = gammas.at(i);
    }
}

void update_alphas (vec& bs_gammas, vec& gammas, vec& alphas,
                    vec& W0H_bs_gammas, vec& W0h_bs_gammas, vec& W0H2_bs_gammas,
                    vec& WH_gammas, vec& Wh_gammas, vec& WH2_gammas,
                    vec& WlongH_alphas, vec& Wlongh_alphas, vec& WlongH2_alphas,
                    vec& log_Pwk, vec& log_Pwk2, uvec& id_H,
                    uvec& which_event, uvec& which_right_event, uvec& which_left, uvec& which_interval,
                    bool& any_event, bool& any_interval,
                    vec& prior_mean_bs_gammas, mat& prior_Tau_bs_gammas, double& tau_bs_gammas,
                    vec& prior_mean_gammas, mat& prior_Tau_gammas,
                    vec& prior_mean_alphas, mat& prior_Tau_alphas,
                    double& denominator_surv, int& it,
                    /////
                    mat& Wlong_H, mat& Wlong_h, mat& Wlong_H2,
                    vec& scale_alphas,
                    mat& acceptance_alphas,
                    mat& res_alphas) {
    for (unsigned int i = 0; i < alphas.n_rows; ++i) {
        vec proposed_alphas = propose_norm(alphas, scale_alphas, i);
        vec proposed_WlongH_alphas = Wlong_H * proposed_alphas;
        vec proposed_Wlongh_alphas(Wlong_h.n_rows);
        if (any_event) {
            proposed_Wlongh_alphas = Wlong_h * proposed_alphas;
        }
        vec proposed_WlongH2_alphas(Wlong_H2.n_rows);
        if (any_interval) {
            proposed_WlongH2_alphas = Wlong_H2 * proposed_alphas;
        }
        double numerator_surv =
            log_density_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                             WH_gammas, Wh_gammas, WH2_gammas,
                             proposed_WlongH_alphas, proposed_Wlongh_alphas, proposed_WlongH2_alphas,
                             log_Pwk, log_Pwk2, id_H,
                             which_event, which_right_event, which_left,
                             any_interval, which_interval) +
            logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas,
                     tau_bs_gammas) +
            logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
            logPrior(proposed_alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
        double log_ratio = numerator_surv - denominator_surv;
        if (is_finite(log_ratio) && exp(log_ratio) > R::runif(0, 1)) {
            alphas = proposed_alphas;
            WlongH_alphas = proposed_WlongH_alphas;
            if (any_event) {
                Wlongh_alphas = proposed_Wlongh_alphas;
            }
            if (any_interval) {
                WlongH2_alphas = proposed_WlongH2_alphas;
            }
            denominator_surv = numerator_surv;
            acceptance_alphas.at(it, i) = 1;
        }
        if (it > 19) {
            scale_alphas.at(i) =
                robbins_monro(scale_alphas.at(i),
                              acceptance_alphas.at(it, i), it);
        }
        // store results
        res_alphas.at(it, i) = alphas.at(i);
    }
}

// [[Rcpp::export]]
List mcmc_cpp (List model_data, List model_info, List initial_values,
               List priors, List control) {
    // outcome vectors and design matrices
    vec Time_right = as<vec>(model_data["Time_right"]);
    vec Time_left = as<vec>(model_data["Time_left"]);
    vec Time_start = as<vec>(model_data["Time_start"]);
    vec delta = as<vec>(model_data["Time_start"]);
    uvec which_event = as<uvec>(model_data["which_event"]) - 1;
    uvec which_right = as<uvec>(model_data["which_right"]) - 1;
    uvec which_right_event = join_cols(which_event, which_right);
    uvec which_left = as<uvec>(model_data["which_left"]) - 1;
    uvec which_interval = as<uvec>(model_data["which_interval"]) - 1;
    mat W0_H = as<mat>(model_data["W0_H"]);
    mat W0_h = as<mat>(model_data["W0_h"]);
    mat W0_H2 = as<mat>(model_data["W0_H2"]);
    mat W_H = as<mat>(model_data["W_H"]);
    mat W_h = as<mat>(model_data["W_h"]);
    mat W_H2 = as<mat>(model_data["W_H2"]);
    mat W_bar = as<mat>(model_data["W_bar"]);
    List Wlong_H_ = as<List>(model_data["Wlong_H"]);
    mat Wlong_H = docall_cbind(Wlong_H_);
    List Wlong_h_ = as<List>(model_data["Wlong_h"]);
    mat Wlong_h = docall_cbind(Wlong_h_);
    List Wlong_H2_ = as<List>(model_data["Wlong_H2"]);
    mat Wlong_H2 = docall_cbind(Wlong_H2_);
    List Wlong_bar_ = as<List>(model_data["Wlong_bar"]);
    field<mat> Wlong_bar = List2Field_mat(Wlong_bar_);
    // other information
    uvec idT = as<uvec>(model_data["idT"]) - 1;
    vec log_Pwk = as<vec>(model_data["log_Pwk"]);
    vec log_Pwk2 = as<vec>(model_data["log_Pwk2"]);
    uvec id_H = create_fast_ind(as<uvec>(model_data["id_H"]));
    uvec id_H2 = create_fast_ind(as<uvec>(model_data["id_H2"]));
    bool any_gammas = as<bool>(model_data["any_gammas"]);
    bool any_event = which_event.n_rows > 0;
    bool any_interval = which_interval.n_rows > 0;
    // initial values
    vec bs_gammas = as<vec>(initial_values["bs_gammas"]);
    vec gammas = as<vec>(initial_values["gammas"]);
    vec alphas = as<vec>(initial_values["alphas"]);
    double tau_bs_gammas = as<double>(initial_values["tau_bs_gammas"]);
    List b_ = as<List>(initial_values["b"]);
    field<mat> b = List2Field_mat(b_);
    mat b_mat = docall_cbind(b_);
    mat D = as<mat>(initial_values["D"]);
    vec sds = sqrt(D.diag());
    mat R = cov2cor(D);
    mat L = chol(R);
    List betas_ = as<List>(initial_values["betas"]);
    field<vec> betas = List2Field_vec(betas_);
    // indexes or other useful things
    uvec upper_part = trimatu_ind(size(R),  1);
    // MCMC settings
    int n_iter = as<int>(control["n_iter"]);
    int n_burnin = as<int>(control["n_burnin"]);
    bool MALA = as<bool>(control["MALA"]);
    // priors
    vec prior_mean_bs_gammas = as<vec>(priors["mean_bs_gammas"]);
    mat prior_Tau_bs_gammas = as<mat>(priors["Tau_bs_gammas"]);
    vec prior_mean_gammas = as<vec>(priors["mean_gammas"]);
    mat prior_Tau_gammas = as<mat>(priors["Tau_gammas"]);
    vec prior_mean_alphas = as<vec>(priors["mean_alphas"]);
    mat prior_Tau_alphas = as<mat>(priors["Tau_alphas"]);
    double post_A_tau_bs_gammas = as<double>(priors["A_tau_bs_gammas"]) +
        0.5 * as<double>(priors["rank_Tau_bs_gammas"]);
    double prior_B_tau_bs_gammas = as<double>(priors["B_tau_bs_gammas"]);
    double prior_D_sds_df = as<double>(priors["prior_D_sds_df"]);
    double prior_D_sds_sigma = as<double>(priors["prior_D_sds_sigma"]);
    double prior_D_L_etaLKJ = as<double>(priors["prior_D_L_etaLKJ"]);
    // store results
    int n_bs_gammas = bs_gammas.n_rows;
    int n_gammas = gammas.n_rows;
    int n_alphas = alphas.n_rows;
    int n_sds = sds.n_rows;
    int n_L = vec(L(upper_part)).n_rows;
    mat res_bs_gammas(n_iter, n_bs_gammas);
    mat acceptance_bs_gammas(n_iter, n_bs_gammas, fill::zeros);
    mat res_gammas(n_iter, n_gammas);
    vec res_W_bar_gammas(n_iter);
    mat acceptance_gammas(n_iter, n_gammas, fill::zeros);
    mat res_alphas(n_iter, n_alphas);
    mat acceptance_alphas(n_iter, n_alphas, fill::zeros);
    vec res_tau_bs_gammas(n_iter);
    mat res_sds(n_iter, n_sds, fill::zeros);
    mat acceptance_sds(n_iter, n_sds, fill::zeros);
    mat res_L(n_iter, n_L, fill::zeros);
    mat acceptance_L(n_iter, n_L, fill::zeros);
    // scales
    vec scale_bs_gammas = create_init_scale(n_bs_gammas);
    vec scale_gammas = create_init_scale(n_gammas);
    vec scale_alphas = create_init_scale(n_alphas);
    vec scale_sds = create_init_scale(n_sds);
    vec scale_L = create_init_scale(n_L);
    // preliminaries
    vec W0H_bs_gammas = W0_H * bs_gammas;
    vec W0h_bs_gammas(W0_h.n_rows);
    vec W0H2_bs_gammas(W0_H2.n_rows);
    if (any_event) {
        W0h_bs_gammas = W0_h * bs_gammas;
    }
    if (any_interval) {
        W0H2_bs_gammas = W0_H2 * bs_gammas;
    }
    vec WH_gammas(W0_H.n_rows);
    vec Wh_gammas(W0_h.n_rows);
    vec WH2_gammas(W0_H2.n_rows);
    if (any_gammas) {
        WH_gammas = W_H * gammas;
    }
    if (any_gammas && any_event) {
        Wh_gammas = W_h * gammas;
    }
    if (any_gammas && any_interval) {
        WH2_gammas = WH2_gammas * gammas;
    }
    vec WlongH_alphas = Wlong_H * alphas;
    vec Wlongh_alphas(W0_h.n_rows);
    if (any_event) {
        Wlongh_alphas = Wlong_h * alphas;
    }
    vec WlongH2_alphas(W0_H2.n_rows);
    if (any_interval) {
        WlongH2_alphas = Wlong_H2 * alphas;
    }
    double denominator_surv =
        log_density_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                         WH_gammas, Wh_gammas, WH2_gammas,
                         WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                         log_Pwk, log_Pwk2, id_H,
                         which_event, which_right_event, which_left,
                         any_interval, which_interval) +
        logPrior(bs_gammas, prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas) +
        logPrior(gammas, prior_mean_gammas, prior_Tau_gammas, 1.0) +
        logPrior(alphas, prior_mean_alphas, prior_Tau_alphas, 1.0);
    for (int it = 0; it < n_iter; ++it) {
        update_bs_gammas(bs_gammas, gammas, alphas,
            W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
            WH_gammas, Wh_gammas, WH2_gammas,
            WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
            log_Pwk, log_Pwk2, id_H,
            which_event, which_right_event, which_left, which_interval,
            any_event, any_interval,
            prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas,
            prior_mean_gammas, prior_Tau_gammas,
            prior_mean_alphas, prior_Tau_alphas,
            denominator_surv, it,
            /////
            W0_H, W0_h, W0_H2, scale_bs_gammas, acceptance_bs_gammas,
            res_bs_gammas);
        ////////////////////////////////////////////////////////////////////////
        double post_B_tau_bs_gammas = prior_B_tau_bs_gammas +
            0.5 * as_scalar(bs_gammas.t() * prior_Tau_bs_gammas * bs_gammas);
        tau_bs_gammas = R::rgamma(post_A_tau_bs_gammas, 1 / post_B_tau_bs_gammas);
        res_tau_bs_gammas.at(it) = tau_bs_gammas;
        ////////////////////////////////////////////////////////////////////////
        if (any_gammas) {
            update_gammas(bs_gammas, gammas, alphas,
                          W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                          WH_gammas, Wh_gammas, WH2_gammas,
                          WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                          log_Pwk, log_Pwk2, id_H,
                          which_event, which_right_event, which_left, which_interval,
                          any_event, any_interval,
                          prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas,
                          prior_mean_gammas, prior_Tau_gammas,
                          prior_mean_alphas, prior_Tau_alphas,
                          denominator_surv, it,
                          /////
                          W_H, W_h, W_H2, scale_gammas, acceptance_gammas,
                          res_gammas);
            res_W_bar_gammas.at(it) = as_scalar(W_bar * gammas);
        }
        ////////////////////////////////////////////////////////////////////////
        update_alphas(bs_gammas, gammas, alphas,
                      W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas,
                      WH_gammas, Wh_gammas, WH2_gammas,
                      WlongH_alphas, Wlongh_alphas, WlongH2_alphas,
                      log_Pwk, log_Pwk2, id_H,
                      which_event, which_right_event, which_left, which_interval,
                      any_event, any_interval,
                      prior_mean_bs_gammas, prior_Tau_bs_gammas, tau_bs_gammas,
                      prior_mean_gammas, prior_Tau_gammas,
                      prior_mean_alphas, prior_Tau_alphas,
                      denominator_surv, it,
                      /////
                      Wlong_H, Wlong_h, Wlong_H2, scale_alphas,
                      acceptance_alphas, res_alphas);
        ////////////////////////////////////////////////////////////////////////
        update_D(L, sds, b_mat, upper_part,
                 prior_D_sds_df, prior_D_sds_sigma, prior_D_L_etaLKJ,
                 it, MALA, res_sds, res_L, scale_sds, scale_L,
                 acceptance_sds, acceptance_L);
    }
    return List::create(
        Named("mcmc") = List::create(
            Named("bs_gammas") = res_bs_gammas.rows(n_burnin, n_iter - 1),
            Named("tau_bs_gammas") = res_tau_bs_gammas.rows(n_burnin, n_iter - 1),
            Named("gammas") = res_gammas.rows(n_burnin, n_iter - 1),
            Named("W_bar_gammas") = res_W_bar_gammas.rows(n_burnin, n_iter - 1),
            Named("alphas") = res_alphas.rows(n_burnin, n_iter - 1),
            Named("sds") = res_sds.rows(n_burnin, n_iter - 1),
            Named("L") = res_L.rows(n_burnin, n_iter - 1)
        ),
        Named("acc_rate") = List::create(
            Named("bs_gammas") = acceptance_bs_gammas.rows(n_burnin, n_iter - 1),
            Named("gammas") = acceptance_gammas.rows(n_burnin, n_iter - 1),
            Named("alphas") = acceptance_alphas.rows(n_burnin, n_iter - 1),
            Named("sds") = acceptance_sds.rows(n_burnin, n_iter - 1),
            Named("L") = acceptance_L.rows(n_burnin, n_iter - 1)
        )
    );
}

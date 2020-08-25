#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

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

field<uvec> List2Field_uvec (const List &uVecs) {
    int n_list = uVecs.size();
    field<uvec> res(n_list);
    for (int i = 0; i < n_list; ++i) {
        res.at(i) = as<uvec>(uVecs[i]);
    }
    return res;
}

uvec create_fast_ind (const uvec &group) {
    unsigned int l = group.n_rows;
    uvec ind = find(group.rows(1, l - 1) != group.rows(0, l - 2));
    unsigned int k = ind.n_rows;
    ind.insert_rows(k, 1);
    ind.at(k) = l - 1;
    return ind;
}

vec group_sum (const vec &x, const uvec &ind) {
    vec cumsum_x = cumsum(x);
    vec out = cumsum_x.elem(ind);
    out.insert_rows(0, 1);
    out = diff(out);
    return out;
}

vec mu_fun (const vec &eta, const std::string &link) {
    uword n = eta.n_rows;
    vec out(n);
    if (link == "identity") {
        out = eta;
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

vec log_long (const field<mat> &y, const field<vec> &eta, const vec &scales,
              const vec &extra_parms, const CharacterVector &families,
              const CharacterVector &links, const field<uvec> &ids,
              const field<uvec> &unq_ids) {
    uword n_outcomes = y.size();
    uvec ns(n_outcomes);
    for (uword i = 0; i < n_outcomes; ++i) {
        ns.at(i) = ids.at(i).n_rows;
    }
    uword n = ns.max();
    vec out(n, fill::zeros);
    for (uword i = 0; i < n_outcomes; ++i) {
        uvec id_i = ids.at(i);
        uvec unq_id_i = unq_ids.at(i);
        mat y_i = y.at(i);
        uword N = y_i.n_rows;
        vec log_contr(N);
        vec mu_i = mu_fun(eta.at(i), std::string(links[i]));
        double scale_i = scales.at(i);
        double extr_prm_i = extra_parms.at(i);
        if (families[i] == "gaussian") {
            log_contr = log_dnorm(y_i, mu_i, scale_i);
        } else if (families[i] == "binomial") {
            uword k = y_i.n_cols;
            if (k == 2) {
                // y_i.col(1) has been set to y_i.col(0) + y_i.col(1)
                // in jm_fit(), i.e., y_i.col(1) is already the number of trials
                // not the number of failures
                log_contr = log_dbinom(y_i.col(0), y_i.col(1), mu_i);
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
                // y_i.col(1) has been set to y_i.col(0) + y_i.col(1)
                // in jm_fit(), i.e., y_i.col(1) is already the number of trials
                // not the number of failures
                log_contr = log_dbbinom(y_i.col(0), y_i.col(1), mu_i, scale_i);
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
        out.rows(unq_id_i) += group_sum(log_contr, id_i);
    }
    return out;
}

// [[Rcpp::export]]
List test_log_long (List model_data, List model_info, List initial_values) {
    List y_ = as<List>(model_data["y"]);
    const field<mat> y = List2Field_mat(y_);
     List eta_ = as<List>(initial_values["eta"]);
    const field<vec> eta = List2Field_vec(eta_);
    const vec scales = exp(as<vec>(initial_values["log_sigmas"]));
    const vec extra_parms(scales.n_rows, fill::zeros);
    CharacterVector families = as<CharacterVector>(model_info["family_names"]);
    CharacterVector links = as<CharacterVector>(model_info["links"]);
    List idL_lp_ = as<List>(model_data["idL_lp"]);
    field<uvec> idL_lp = List2Field_uvec(idL_lp_);
    for (uword i = 0; i < idL_lp.size(); ++i) {
        idL_lp.at(i) = create_fast_ind(idL_lp.at(i));
    }
    List unq_idL_ = as<List>(model_data["unq_idL"]);
    field<uvec> unq_idL = List2Field_uvec(unq_idL_);
    for (uword i = 0; i < unq_idL.size(); ++i) {
        unq_idL.at(i) = unq_idL.at(i) - 1;
    }
    vec out = log_long(y, eta, scales, extra_parms, families, links, idL_lp,
                       unq_idL);
    return List::create(
        Named("y") = y,
        Named("scales") = scales,
        Named("extra_parms") = extra_parms,
        Named("idL_lp") = idL_lp,
        Named("families") = families,
        Named("links") = links,
        Named("unq_idL") = unq_idL,
        Named("log_Lik") = out
    );
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

// [[Rcpp::export]]
List test (List model_data, List model_info, List initial_values) {
    List X_H_ = as<List>(model_data["X_H"]);
    const field<mat> X_H = List2Field_mat(X_H_);
    List X_h_ = as<List>(model_data["X_h"]);
    const field<mat> X_h = List2Field_mat(X_H_);
    List X_H2_ = as<List>(model_data["X_H2"]);
    const field<mat> X_H2 = List2Field_mat(X_H2_);
    List Z_H_ = as<List>(model_data["Z_H"]);
    const field<mat> Z_H = List2Field_mat(Z_H_);
    List Z_h_ = as<List>(model_data["Z_h"]);
    const field<mat> Z_h = List2Field_mat(X_H_);
    List Z_H2_ = as<List>(model_data["Z_H2"]);
    const field<mat> Z_H2 = List2Field_mat(Z_H2_);
    List b_ = as<List>(initial_values["b"]);
    field<mat> b = List2Field_mat(b_);
    List betas_ = as<List>(initial_values["betas"]);
    field<vec> betas = List2Field_vec(betas_);
    //uvec id_H = create_fast_ind(as<uvec>(model_data["id_H"]));
    uvec id_H = as<uvec>(model_data["id_H"]) - 1;
    field<mat> eta_H = linpred_surv(X_H, betas, Z_H, b, id_H);
    return List::create(
        Named("X_H") = X_H,
        Named("b") = b,
        Named("betas") = betas,
        Named("eta_H") = eta_H
    );
}

field<mat> create_Wlong (const field<mat> &eta, const field<uvec> &FunForms,
                         const field<mat> &U, const field<uvec> &ind) {
    uword n_outcomes = eta.n_elem;
    field<mat> out(n_outcomes);
    for (uword i = 0; i < n_outcomes; ++i) {
        mat eta_i = eta.at(i);
        uvec FF_i = FunForms.at(i);
        mat U_i = U.at(i);
        uvec ind_i = ind.at(i);
        mat Wlong_i(eta_i.n_rows, U_i.n_cols, fill::ones);
        Wlong_i.cols(FF_i) %= eta_i.cols(ind_i);
        out.at(i) = U_i % Wlong_i;
    }
    return out;
}



#ifndef JMBAYES2FE_H
#define JMBAYES2FE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

/*?? To update on jm()
ind_FE_nHC <- which(!seq_len(sum(nfes)) %in% ind_FE_HC)  # new
indL_FE_nHC <- mapply2(function (x, ind) seq_along(x)[-ind], ind_FE, x_in_z_base) # new
has_tilde_betas <- sapply(ind_FE_nHC, length) > 0
indL_FE_nHC[] <- lapply(ind_FE_nHC, function (x) if (length(x)) x else 0L) # new because changed the name of the var, the remaining is the same
OT: ...[] what is the benefit of using the brackets? is it fast because because updates only the content and keeps the old attributes?
*/ 
  
/*?? To update on jm()
change q_dot to nfes_HC 
*/

/*?? To update on C functions:
field<vec> vec2field (const vec &betas, const field<uvec> &ind_FE) {
  uword n = ind_FE.n_elem;
  field<vec> out(n);
  for (uword i = 0; i < n; i++) {
    out.at(i) = betas.rows(ind_FE.at(i) - 1); # <- The change is here, with the "-1". The counter in ind_FE starts from 1, and not 0. I cannot do ind_FE-1, because ind_FE is a field.
  }
  return out;
}
*/

void update_betas (field<vec> &betas, // it-th sampled fixed effects
                   mat &res_betas, // all sampled fixed effects
                   const uword &it, // current iteration
                   const field<vec> &prior_mean_FE, const field<mat> &prior_Tau_FE, // prior
                   const mat &b_mat; // it-th sampled RE
                   const mat &L, // RE corr matrix factorization (lower matrix)
                   const vec &sds, // RE SDs
                   const mat &X_dot, // X_dot matrix
                   const field<uvec> &ind_FE // indices for the FE in res_betas[it,] belonging to the field betas. E.g., {{1,2,3}, {4, 5}, {6}}
                   const uvec &ind_FE_HC, // indices for the FE present in the HC (cols in res_betas). Counts from 1.
                   const uvec &ind_FE_nHC, // indices for the FE NOT present in the HC (cols in res_betas) Counts from 1.
                   const field<uvec> &indL_FE_nHC, // indices for the FE NOT present in the HC per outcome (rows is each field element. the count is reset for each outcome). Counts from 1.
                   const vec &id_patt, // vector with the ids' outcome missing pattern
                   const field<uvec> &ind_RE_patt, // indices for the RE present in each outcome missing pattern (cols in D). Counts from 1.
                   const field<uvec> &ind_FE_patt, // indices for the FE (in HC) present in each outcome missing pattern (cols in X_dot). Counts from 1.
                   mat &acceptance_betas,
                   vec &scale_betas,
                   const uvec &has_tilde_betas,
                   const field<mat> &X, 
                   const field<mat> &Z,
                   const field<mat> &b,
                   const field<uvec> &id,
                   const field<mat> &y,
                   const vec &sigmas,
                   const vec &extra_parms,
                   const CharacterVector &families,
                   const CharacterVector &links,
                   const field<uvec> &ids,
                   const field<uvec> &unq_ids,
                   const field<mat> &X_H, const field<mat> &X_h, const field<mat> &X_H2,
                   const field<mat> &Z_H, const field<mat> &Z_h, const field<mat> &Z_H2,
                   const field<mat> &U_H, const field<mat> &U_h, const field<mat> &U_H2,
                   const mat &Wlong_bar,
                   const uvec &id_H_, const uvec &id_h,
                   const field<uvec> &FunForms,
                   const field<uvec> &FunForms_ind,
                   const vec &alphas,
                   mat &Wlong_h, mat &Wlong_H, mat &Wlong_H2,
                   vec &Wlongh_alphas,
                   const bool &any_event, const bool &any_interval,
                   const vec &W0H_bs_gammas, const vec &W0h_bs_gammas, const vec &W0H2_bs_gammas,
                   const vec &WH_gammas, const vec &Wh_gammas, const vec &WH2_gammas,
                   const vec &log_Pwk, const vec &log_Pwk2,
                   const uvec &indFast_H, const uvec &indFast_h,
                   const uvec &which_event, const uvec &which_right_event, const uvec &which_left,
                   const uvec &which_interval,
                   const field<uvec> prop_Sigma_FE_nHC, // vcov for the proposal
                   vec &logLik_long,
                   vec &logLik_surv,
                   vec &logLik_fe_nHC, // a vec(n_outcomes) where vec.rows(j) is the logLik_prior for outcome-j at iteration it-1  
                   const uword &RMu_it_thld, // robins monro univariate iterations threshold to update the scale. Updates scale from it > RMu_it_thld
                   const double &RMu_tacce // robins monro univariate target acceptance
                  )
   
  // For the HC-FE we first update the row on the res_betas and then update the field betas 
  // For the nHC-FE outside HC we update the field betas per outcome, and after updating all outcomed we update the row in res_betas

  // FE in HC - Gibss sampling
  
  //?? The variables below could also be function's input parameters. They not change between iterations.
  uword n = b_mat.n_rows; // number of unique subjects
  uword n_y = betas.n_elem; // number of longitudinal outcomes
  uword re_count = L.n_cols; // number of random effects
  uword patt_count = id_patt.max(); // number of unique outcome-missing patterns          
  // below: it assumes that the precision matrices within and between outcomes are diagonal. If not, we need to implement a function to obtain the conditional prior distributions.
  vec prior_mean_HC = docall_rbindF(prior_mean_FE).rows(ind_FE_HC - 1); // prior mean for all HC-FE
  mat prior_Tau_HC  = bdiagF(prior_Tau_FE).submat(ind_FE_HC - 1, ind_FE_HC - 1); // prior precision for all HC-FE
  uword p_hc = ind_FE_HC.n_elem; // number of HC-FE
  
  mat sum_JXDXJ(p_hc, p_hc, fill::zeros); // sum for the posterior parameters
  vec sum_JXDu(p_hc, fill::zeros); // sum for the posterior parameters
  mat L_D = L.each_row() % sds.t(); // RE vcov matrix Choesky factorization (upper)
  //?? the L matrix is a lower or upper triangular? If upper update the L_... variables below
  vec betas_HC(p_hc) = docall_rbindF(betas).rows(ind_FE_HC - 1);
  /*?? I am declaring the variables outside the loop, because it seems to reduce the computation time. 
  Seen here: https://gallery.rcpp.org/articles/dmvnorm_arma/. I am doing the same for the loop
  in FE-nHC. Do you have any thoughts on this? */
  
  uword patt_i; // outcome missing pattern for the i-th patient
  mat L_patt_inv; // inv cholesky factorization vcov matrix for the i-th (<patt_count) missing pattern
  field<mat> D_inv(patt_count); // all unique vcov_inv matrices accross the missing outcome patterns
  mat X_dot_i; // X_dot for the i-th patient
  vec u_i; // all RE for the i-th patient (including the RE for the missing outcomes. they'll be cancelled out by the zero-row)
  mat D_inv_i; // inverse vcov matrix for the i-th patient
  mat XD_i; // aux var for the posterior parameters for the i-th patient
  mat XDX_i; // aux var for the posterior parameters for the i-th patient
  uvec ind_FE_i; // ind of the HC-FE for the i-th patient (patients reporting different outcomes will have different ind)
  
  vec logLik_long_prop; // all longitudinal outcomes contribution for the logLik (given the proposed nHC-FE). each k-th row reports to the k-th id
  vec logLik_surv_prop; // survival contribution for the logLik (given the proposed nHC-FE). each k-th row reports to the k-th id
  mat Wlong_H_prop;
  vec WlongH_alphas_prop; 
  mat Wlong_h_prop(Wlong_h.n_rows, Wlong_h.n_cols);
  vec Wlongh_alphas_prop(Wlongh_alphas.n_rows);
  mat Wlong_H2_prop(Wlong_H2.n_rows, Wlong_H2.n_cols);
  vec WlongH2_alphas_prop(WlongH2_alphas.n_rows);
  
  for (uword i = 0; i < n; ++i) { // i-th patient
    
    patt_i = id_patt.at(i); // id missing outcome pattern // uword
    ind_FE_i = ind_FE_patt.at(patt_i) - 1; // uvec
    
    if (i < patt_count) { // obtain all unique vcov_inv matrices required for the sums in the posterior parameters
      
      L_patt_inv = inv( trimatl( chol_update(L_D, ind_RE_patt.at(i) - 1) ) ); // mat
      D_inv.at(i) =  L_patt_inv.t() * L_patt_inv; // mat //?? do you know a better way to do this caculation?
      
    }

    X_dot_i = X_dot.submat(ind_RE_patt.at(patt_i) - 1 + i*re_count, 
                           ind_FE_i); // mat
    
    u_i = b_mat.row(i).t() + X_dot_i * betas_HC.elem(ind_FE_i); // vec
    D_inv_i = D_inv.at(patt_i); // mat 
    
    XD_i = X_dot_i.t() * D_inv_i; // mat
    XDX_i = XD_i * X_dot_i; // mat
    
    sum_JXDu  += add_zero_colrows(XD_i, p_hc, XD_i.n_cols, ind_FE_i, regspace<uvec>(0,  XD_i.n_cols - 1)) * u_i; // mat
    sum_JXDXJ += add_zero_colrows(XDX_i, p_hc, p_hc, ind_FE_i, ind_FE_i); // mat
    
  }
  
  mat Tau_1 = inv(prior_Tau_HC + sum_JXDXJ); //?? Can I write the Tau_1 in terms of an L matrix? to avoid the use of chol(Tau_1) later <---------------------------------
  vec mean_1 = Tau_1 * (prior_Tau_HC * prior_mean_HC + sum_JXDu);
  mat U_1 = chol(Tau_1);
  
  res_betas.submat(it, ind_FE_HC - 1) = (propose_mvnorm_vec(1, U_1, 1) + mean_1).t(); // vec.t()
  /* below: I am doing this here, because I need to the save the update the sampled FE_HC 
  on the betas field, and I need to add the current FE_nHC to avoid add zeros instead. 
  Given I need the current FE values for the MH. I later update the FE_nHC on this 
  same row with the M-H sample. By doing so we avoid a create new ind_variables. */
  res_betas.submat(it, ind_FE_nHC - 1) = docall_rbindF(betas).rows(ind_FE_nHC - 1).t(); // vec.t()
  betas = vec2field(res_betas.row(it), ind_FE); // field<vec>
  
  //?? I believe I need to update the logLik_long, logik_surv, and Wlongh... after each iteration of the Gibbs sampler. Right?
  
  // logLik_long
  eta_long = linpred_mixed(X, betas, Z, b, id); // field<vec>
  logLik_long_prop = log_long(y, eta_long, sigmas, extra_parms,
                              families, links, ids, unq_ids, n); // vec
  /* there is no need to estimate the logLik_long_prop for all outcomes, the logLik its proportional to the contibution of outcome j. 
  we could remove the raining outcomes from the posterior. But we are using all outcomes 
  because if we accept we will need to update the (full) logLik_long variable */
  
  // logLik_surv
  Wlong_H_prop = calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, betas, b, 
                                 id_H_, FunForms, FunForms_ind); // mat
  WlongH_alphas_prop = Wlong_H_prop * alphas; // vec
  
  if (any_event) {
    Wlong_h_prop = calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, betas, b,
                                   id_h, FunForms, FunForms_ind); // mat
    Wlongh_alphas_prop = Wlong_h_prop * alphas; // vec
  }
  
  if (any_interval) {
    Wlong_H2_prop = calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, betas, b, 
                                    id_H_, FunForms, FunForms_ind); // mat
    WlongH2_alphas_prop = Wlong_H2_prop * alphas; // vec
  }
  
  
  logLik_surv_prop = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas, 
                              WH_gammas, Wh_gammas, WH2_gammas, 
                              WlongH_alphas_prop, Wlongh_alphas_prop,
                              WlongH2_alphas_prop, log_Pwk, log_Pwk2, 
                              indFast_H, indFast_h, which_event, which_right_event, 
                              which_left, any_interval, which_interval); // vec
  
  Wlong_H = Wlong_H_prop;
  WlongH_alphas = WlongH_alphas_prop;
  
  if (any_event) {
    Wlong_h = Wlong_h_prop;
    Wlongh_alphas = Wlongh_alphas_prop;
  }
  
  if (any_interval) {
    Wlong_H2 = Wlong_H2_prop;
    WlongH2_alphas = WlongH2_alphas_prop;
  }
  
  //////////////////////////////////////////////////////////////////////////////
  // FE not in HC - Metropolis-Hastings sampling
  
  uvec ind_j; // ind nHC-FE in the j-th outcome. Counts from 0 
  mat prop_U_j; // Chol-fact (upper) vcov matrix proposal distribution for the nHC-FE in the j-th outcome
  vec prior_mean_j; // prior mean for the nHC-FE in the j-th outcome
  mat prior_U_j; //  Chol-fact (upper) vcov matrix for the nHC-FE in the j-th outcome
  field<vec> eta_long; // eta for all longitudinal outcomes
  double denominator; // log acceptance ratio denominator
  double numerator; // log acceptance ratio numerator 
  double log_ratio; // log acceptance ratio
  field<uvec> betas_prop; // all FE with the proposed nHC-FE for the j-th outcome
  double logLik_fe_nHC_prop; // prior contribution from the j-th outcome for the logLik (given the proposed nHC-FE)
  
  for (uword j = 0; j < n_y; ++i) { // j-th outcome
    
    if (has_tilde_betas.rows(it)) { // skip the outcomes that don't have nHC-FE
      
      ind_j = indL_FE_nHC.at(j) - 1; // uvec
      //?? the three lines below could be outside the function. They do not change accross iterations
      prop_U_j = chol(prop_Sigma_FE_nHC.at(j).submat(ind_j, ind_j)); // mat
      prior_mean_j = prior_mean_FE.at(j).rows(ind_j);  // vec
      // below: log_dmvnrm_chol() expects a chol(Sigma) and not a chol(Tau). It assumes that prior_Tau_FE.at(j) is a diagonal matrix.
      prior_U_j = chol( inv(diagmat( prior_Tau_FE.at(j).submat(ind_j, ind_j) )) ); // mat 
      /*?? By using chol(.) we avoid uncessary chol(.) calculations. If we were 
      using log_dmvnrm() instead, it will do a chol(Sigma) for each iteration 
      and outcome. By using log_dmvnrm_chol(), we require only one chol() per outcome. 
      But this is only true if we create these variables outside this function.*/
      
      if(it = 0) { //?? Greg seems to do this step outside. Check this with him.
        logLik_fe_nHC.rows(j) = log_dmvnrm_chol(betas.at(j).rows(ind_j) - prior_mean_j,
                                              prior_U_j); // double
        }

      // denominator
      denominator = sum(logLik_long) + sum(logLik_surv) + logLik_fe_nHC.at(j); // double
 
      betas_prop = betas; // field<vec>
      beta_prop.at(j).rows(ind_j) = propose_mvnorm_vec(1, pro_U_j, scale_betas.rows(j)) + betas.at(j).rows(ind_j); // vec
      
      // logLik_prior
      logLik_fe_nHC_prop = log_dmvnrm_chol(betas_prop.at(j).rows(ind_j) - prior_mean_j,
                                           prior_U_j); // double

      // logLik_long
      eta_long = linpred_mixed(X, betas_prop, Z, b, id); // field<vec>
      logLik_long_prop = log_long(y, eta_long, sigmas, extra_parms,
                                  families, links, ids, unq_ids, n); // vec
      /* there is no need to estimate the logLik_long_prop for all outcomes, the 
      logLik its proportional to the contibution of outcome j. 
      We could remove the raining outcomes from the posterior. But we are using 
      all outcomes because if we accept we will need to update the (full) 
      logLik_long variable */
      
      // logLik_surv
      Wlong_H_prop = calculate_Wlong(X_H, Z_H, U_H, Wlong_bar, betas_prop, b, 
                                     id_H_, FunForms, FunForms_ind); // mat
      WlongH_alphas_prop = Wlong_H_prop * alphas; // vec
      
      if (any_event) {
        Wlong_h_prop = calculate_Wlong(X_h, Z_h, U_h, Wlong_bar, betas_prop, b,
                                       id_h, FunForms, FunForms_ind); // mat
        Wlongh_alphas_prop = Wlong_h_prop * alphas; // vec
      }
      
      if (any_interval) {
        Wlong_H2_prop = calculate_Wlong(X_H2, Z_H2, U_H2, Wlong_bar, betas_prop, b, 
                                        id_H_, FunForms, FunForms_ind); // mat
        WlongH2_alphas_prop = Wlong_H2_prop * alphas; // vec
      }
      
      logLik_surv_prop = log_surv(W0H_bs_gammas, W0h_bs_gammas, W0H2_bs_gammas, 
                                  WH_gammas, Wh_gammas, WH2_gammas, 
                                  WlongH_alphas_prop, Wlongh_alphas_prop,
                                  WlongH2_alphas_prop, log_Pwk, log_Pwk2, 
                                  indFast_H, indFast_h, which_event, which_right_event, 
                                  which_left, any_interval, which_interval); // vec
      
      // numerator
      numerator = sum(logLik_long_prop) + sum(logLik_surv_prop) + logLik_fe_nHC_prop; // double
      
      /*?? If I'm understanding Greg's code correctly, he samples all patients RE at the
      same time and then updates the RE for those patients who are accepted. 
      He does it because in his case the posterior is proportional to one subject at a time.
      That's why he can sample all at the same time, without the need to conditional on the previous update.
      Did i understand it right? */ 
      
      log_ratio =  numerator - denominator; // double
      
      if(std::isfinite(log_ratio) && std::exp(log_ratio) > R::runif(0.0, 1.0)) {
        
        betas.at(j).rows(ind_j) = betas_prop.at(j).rows(ind_j); // vec
        acceptance_betas.at(it, j) = 1;

        logLik_long = logLik_long_prop; // vec
        logLik_surv = logLik_surv_prop; // vec
        logLik_fe_nHC.rows(j) = logLik_fe_nHC_prop; // double
        
        //?? I would to cross-check the updates below
        Wlong_H = Wlong_H_prop;
        WlongH_alphas = WlongH_alphas_prop;
        
        if (any_event) {
          Wlong_h = Wlong_h_prop;
          Wlongh_alphas = Wlongh_alphas_prop;
        }
        
        if (any_interval) {
          Wlong_H2 = Wlong_H2_prop;
          WlongH2_alphas = WlongH2_alphas_prop;
        }
      }
      
      if(it > RMu_it_thld) {
        scale_betas.rows(j) = robbins_monro(scale_betas.rows(j), acceptance_betas.at(it, j), it, RMu_tacce); // double
      }
    }
  }
  
  res_betas.submat(it, ind_FE_nHC - 1) = docall_rbindF(betas).rows(ind_FE_nHC - 1).t(); // vec.t()
  
)

#endif

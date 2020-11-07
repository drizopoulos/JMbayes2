#ifndef JMBAYES2FE_H
#define JMBAYES2FE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// To update on jm()
  // ind_FE_nHC <- which(!seq_len(sum(nfes)) %in% ind_FE_HC)  # new
  // indL_FE_nHC <- mapply2(function (x, ind) seq_along(x)[-ind], ind_FE, x_in_z_base) # new
  // has_tilde_betas <- sapply(ind_FE_nHC, length) > 0
  // indL_FE_nHC[] <- lapply(ind_FE_nHC, function (x) if (length(x)) x else 0L) # new because changed the name of the var, the remaining is the same

// To update on jm()
  // change q_dot to nfes_HC 

// To update on C functions:
  // field<vec> vec2field (const vec &betas, const field<uvec> &ind_FE) {
  //   uword n = ind_FE.n_elem;
  //   field<vec> out(n);
  //   for (uword i = 0; i < n; i++) {
  //     out.at(i) = betas.rows(ind_FE.at(i) - 1); # <- The change is here, with the "-1". The counter in ind_FE starts from 1, and not 0. I cannot do ind_FE-1, because ind_FE is a field.
  //   }
  //   return out;
  // }

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
                   const uword &RMu_it_thld, // robins monro univariate iterations threshold to update the scale. Updates scale from it > RMu_it_thld
                   const double &RMu_tacce // robins monro univariate target acceptance
                  )
   
   
  // The variables below could also be function's input parameters. They not change between iterations.
  uword n = b_mat.n_rows; // number of unique subjects
  uword n_y = betas.n_elem; // number of longitudinal outcomes
  uword re_count = L.n_cols; // number of random effects
  uword patt_count = id_patt.max(); // number of unique outcome-missing patterns          
  
  vec joint_prior_mean_FE2 = docall_rbindF(prior_mean_FE); // bind the mean-priors from all outcomes in one vector
  mat joint_prior_Tau_FE2  = bdiagF(prior_Tau_FE);  // bind the Tau-priors from all outcomes in one block-matrix
  
  vec prior_mean_FE_HC = joint_prior_mean_FE.rows(ind_FE_HC - 1); // prior for FE in HC 
  mat prior_Tau_FE_HC  = joint_prior_Tau_FE.submat(ind_FE_HC - 1, ind_FE_HC - 1);
  
  vec prior_mean_FE_nHC = joint_prior_mean_FE.rows(ind_FE_nHC - 1); // prior for FE NOT in HC 
  mat prior_Tau_FE_nHC  = joint_prior_Tau_FE.submat(ind_FE_nHC - 1, ind_FE_nHC - 1);
  
  // For the FE in HC we first update the full row on the res_betas and then update the field betas. 
  // For the FE outside HC we update the field betas per outcome, and then update the full row in res_betas

  // FE in HC - Gibss sampling
  uword p_hc = ind_FE_HC.n_elem; // number of FE in the HC
  mat sum_JXDXJ(p_hc, p_hc, fill::zeros); // sum required for posterior parameters
  vec sum_JXDu(p_hc, fill::zeros); // sum required for posterior parameters
  
  mat L_D = L.each_row() % sds.t(); // RE vcov matrix factorization
  field<mat> D_inv(patt_count); // all unique vcov_inv matrices accross the missing outcome patterns

  vec betas_HC(p_hc) = docall_rbindF(betas).rows(ind_FE_HC - 1)

  uword patt_i;
  mat L_patt_inv
  mat X_dot_i;
  vec u_i;
  mat D_inv_i;
  mat XD_i;
  mat XDX_i;
  
  for (uword i = 0; i < n; ++i) { // obtain the sums required for the posterior parameters
    
    patt_i = id_patt.at(i); // id missing outcome pattern
    
    if (i < patt_count) { // obtain all unique vcov_inv matrices required for the sums in the posterior parameters
      
      L_patt_inv = inv( trimatl( chol_update(L_D, ind_RE_patt.at(i) - 1) ) );
      D_inv.at(i) =  L_patt_inv.t() * L_patt_inv;
      
    }

    X_dot_i = X_dot.submat(ind_RE_patt.at(patt_i) - 1 + i*re_count, 
                           ind_FE_patt.at(patt_i) - 1); 
    
    u_i = b_mat.row(i).t() + X_dot_i * betas_HC.elem(ind_FE_patt.at(patt_i) - 1);
    D_inv_i = D_inv.at(patt_i);  
    
    XD_i = X_dot_i.t() * D_inv_i;
    XDX_i = XD_i * X_dot_i;
    
    sum_JXDu  += add_zero_colrows(XD_i, p_hc, XD_i.n_cols, ind_FE_patt.at(patt_i) - 1, regspace<uvec>(0,  XD_i.n_cols - 1)) * u_i; // add zero-rows
    sum_JXDXJ += add_zero_colrows(XDX_i, p_hc, p_hc, ind_FE_patt.at(patt_i) - 1, ind_FE_patt.at(patt_i) - 1); // add zero-rows and zero-columns
    
  }
  
  vec mean_1; // posterior mean vector
  mat Tau_1; // posterior vcov matrix
  
  Tau_1 = (prior_Tau_FE_HC + sum_JXDXJ).i(); //! Can I write the Tau_1 in terms of an L matrix? to avoid the use of chol(Tau_1) later <---------------------------------
  mean_1 = Tau_1 * (prior_Tau_FE_HC * prior_mean_FE_HC + sum_JXDu);
  cube U_1(p_hc, p_hc, 1); U_1.slice(0) = chol(Tau_1); // propose_mvnorm_mat() expects a cube
  
  res_betas.submat(it, ind_FE_HC - 1) = (propose_mvnorm_mat(1, U_1, {1}) + mean_1).t();
  
  // we update the FE_HC in betas after the M-H sampling
  
  // FE not in HC - Metropolis-Hastings sampling
  
  // (TEMPORARY) additional function inputs
  const field<uvec> &L_prop_betas
  
  // (TEMPORARY) other variables
  vec proposed_betas;
  double numerator;
  double denominator;
  
  for (uword i = 0; i < n_y; ++i) {
    
    if (has_tilde_betas.rows(it)) { // some outcomes may not have FE out of the HC
      
      if(it = 0) { denominator = target_log_dist(betas.at(i).rows(indL_FE_nHC.at(i) - 1)) }
      
      proposed_betas = MNORM(1, L_prop_betas, scale) + betas.at(i).rows(indL_FE_nHC.at(i) - 1); // UPDATE <-------------------------
      
      numerator = TARGET_LOG_DIST(proposed_betas); // UPDATE <-------------------------
      
      log_ratio =  numerator - denominator;
      
      if(std::isfinite(log_ratio) && std::exp(log_ratio) > R::runif(0.0, 1.0)) {
        
        betas.at(i).rows(indL_FE_nHC.at(i) - 1) = proposed_betas;
        acceptance_betas.at(it, i) = 1;
        denominator = numerator;
        
      }
      
      if(it > RMu_it_thld) {
        scale_betas.rows(i) = robbins_monro(scale_betas.rows(i), acceptance_betas.at(it, i), it, RMu_tacce);
      }
    }
  }
  
  res_betas.submat(it, ind_FE_nHC - 1) = docall_rbindF(betas).rows(ind_FE_nHC - 1);
  
  betas = vec2field(res_betas.row(it), ind_FE);
  
)

#endif

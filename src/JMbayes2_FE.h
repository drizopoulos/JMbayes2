#ifndef JMBAYES2FE_H
#define JMBAYES2FE_H

#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

void update_betas (field<vec> &betas, // it-th sampled fixed effects
                   mat &res_betas, // all sampled fixed effects
                   const uword &it, // current iteration
                   const field<vec> &prior_mean_FE, const field<mat> &prior_Tau_FE, // prior
                   const mat &b_mat; // it-th sampled RE
                   const mat &L, // RE corr matrix factorization (lower matrix)
                   const vec &sds, // RE SDs
                   const mat &X_dot, // X_dot matrix
                   const field<uvec> &ind_FE, // indices for the FE present in each outcome (to save the sampled FE into a field)
                   const uvec &ind_FE_HC, // indices for the FE present in the HC (cols in res_betas)
                   const uvec &ind_FE_nHC, // indices for the FE NOT present in the HC (cols in res_betas)
                   const vec &id_patt, // vector with the ids' outcome missing pattern
                   const field<uvec> &ind_RE_patt, // indices for the RE present in each outcome missing pattern (cols in D)
                   const field<uvec> &ind_FE_patt // indices for the FE (in HC) present in each outcome missing pattern (cols in X_dot)
                  )
   
   
  // The variables below could also be function's input parameters. They not change between iterations.
  uword n = b_mat.n_rows; // number of unique subjects
  uword re_count = L.n_cols; // number of random effects
  uword patt_count = id_patt.max(); // number of unique outcome-missing patterns          
  
  vec prior_mean_FE2 = docall_rbindF(prior_mean_FE); // bind the mean-priors from all outcomes in one vector
  mat prior_Tau_FE2  = bdiagF(prior_Tau_FE);  // bind the Tau-priors from all outcomes in one block-matrix
  
  vec prior_mean_FE_HC = prior_mean_FE2.rows(ind_FE_HC - 1); // prior for FE in HC 
  mat prior_Tau_FE_HC  = prior_Tau_FE2.submat(ind_FE_HC - 1, ind_FE_HC - 1);
  
  vec prior_mean_FE_nHC = prior_mean_FE2.rows(ind_FE_nHC - 1); // prior for FE NOT in HC 
  mat prior_Tau_FE_nHC  = prior_Tau_FE2.submat(ind_FE_nHC - 1, ind_FE_nHC - 1);
  
  // For the FE in HC we update the full row on the res_betas and then update the field betas. 
  // For the FE outside HC we do the opposite, i.e., we update per outcome in the field betas, and then update the full row in res_betas

  // FE in HC - Gibss sampling
  uword p_hc = ind_FE_HC.n_elem; // number of FE in the HC
  mat sum_JXDXJ(p_hc, p_hc, fill::zeros); // sum required for posterior parameters
  vec sum_JXDu(p_hc, fill::zeros); // sum required for posterior parameters
  
  mat L_D = L.each_row() % sds.t(); // RE vcov matrix factorization
  field<mat> D_inv(patt_count); // all unique vcov_inv matrices accross the missing outcome patterns

  if(it = 0) { 
    vec current_FE_HC(p_hc) = docall_rbindF(betas).rows(ind_FE_nHC - 1)
  } else { 
    vec current_FE_HC(p_hc) = res_betas.submat(it, ind_FE_nHC - 1).t(); // skips the docall_rbindF() step
  }
  
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
      
      L_patt_inv = inv( trimatl( chol_update(L_D, ind_RE_patt.at(i)) ) );
      D_inv.at(i) =  L_patt_inv.t() * L_patt_inv;
      
    }

    X_dot_i = X_dot.submat(ind_RE_patt.at(patt_i) - 1 + i*re_count, 
                           ind_FE_patt.at(patt_i) - 1); 
    
    u_i = b_mat.row(i).t() + X_dot_i * current_FE_HC.elem(ind_FE_patt.at(patt_i) - 1);
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
  
  // FE not in HC - Metropolis-Hastings sampling
  
  //////////////////////////////////////////////////////////////////////////////
  // Draft of the required code
  
  // additional function inputs
  mat &acceptance_betas
  const field<uvec> &indF_FE_nHC // indices for the FE NOT present in the HC per outcome
  vec &scale_betas
  const field<uvec> &vcov_prop_betas
  const uword &RMu_it_thld // robins monro univariate iterations threshold to update the scale
  const double &RMu_acce // robins monro univariate target acceptance
  const field<uvec> &ind_FE // needs to be something like {{1,2,3}, {1, 2}, {1}}. It seems to be generated on jm().
  // other variables
  
  uword n_y = ; // number of outcomes with FE not in HC
  vec proposed_betas;
 
  double numerator;
  double denominator;
  
  for (uword i = 0; i < n_y; ++i) { // this is not correct, I need to account for the tfact that some outcome may not need this step
    
    if(it = 0) { denominator = target_log_dist(betas.at(i).rows(indF_FE_nHC.at(i) - 1)) }
    
    proposed_betas = rmvnorm(1, L, scale) + betas.at(i).rows(indF_FE_nHC.at(i) - 1); // use the correct function
    
    numerator = target_log_dist(proposed_betas); // use the correct function
      
    log_ratio =  numerator - denominator;
    
    if(log_ratio > log(runif(1))) {
      betas.at(i).rows(indF_FE_nHC.at(i) - 1) = proposed_betas;
      acceptance_betas.at(i, it) = 1; // check if the way toa cces acceptance betas is correct 
      denominator = numerator;
    }
    
    if(it > RMu_it_thld) {
      scale_betas.at(i) = robbins_monro_univ(scale = scale_betas.at(i), // use the corrrect function
                                           acceptance_it = acceptance_betas.at(i, it),
                                           it = it, target_acceptance = RMu_acce)
    }
    
  }
  
  res_betas.submat(it, ind_FE_nHC - 1) = docall_rbindF(betas).rows(ind_FE_nHC - 1);
  
  //////////////////////////////////////////////////////////////////////////////
  
  // the FE in HC were not update yet inside the field betas. We add the betas and 'again' the betas_tilde. Doing this here makes the code 'lighter' 
  betas = vec2field(res_betas.row(it), ind_FE);
  
)

#endif

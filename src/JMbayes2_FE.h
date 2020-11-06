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
                   const uvec &ind_FE_in_hc, // indices for the FE present in the HC (cols in res_betas)
                   const uvec &ind_FE_notin_hc, // indices for the FE NOT present in the HC (cols in res_betas)
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
  
  vec prior_mean_FE_in_hc = prior_mean_FE2.elem(ind_FE_in_hc - 1); // prior for FE in HC 
  mat prior_Tau_FE_in_hc  = prior_Tau_FE2.submat(ind_FE_in_hc - 1, ind_FE_in_hc - 1);
  
  vec prior_mean_FE_notin_hc = prior_mean_FE2.elem(ind_FE_notin_hc - 1); // prior for FE NOT in HC 
  mat prior_Tau_FE_notin_hc  = prior_Tau_FE2.submat(ind_FE_notin_hc - 1, ind_FE_notin_hc - 1);
  

  // FE in HC - Gibss sampling
  uword p_hc = ind_FE_in_hc.n_elem; // number of FE in the HC
  mat J(p_hc, p_hc, fill::eye); // identity matrix, required in JXDXJ_i and JXDu_i
  mat sum_JXDXJ(p_hc, p_hc, fill::zeros); // sum required for posterior parameters
  vec sum_JXDu(p_hc, fill::zeros); // sum required for posterior parameters
  
  mat L_D = L.each_row() % sds.t(); // RE vcov matrix factorization
  field<mat> D_inv(patt_count, 1); // all unique vcov_inv matrices accross the missing outcome patterns

  if(it = 0) { 
    vec current_FE_in_hc(p_hc) = docall_rbindF(betas).rows(ind_FE_notin_hc - 1)
  } else { 
    vec current_FE_in_hc(p_hc) = res_betas.submat(it, ind_FE_notin_hc - 1).t(); // skips the docall_rbindF() step
  }
  
  for (uword i = 0; i < n; ++i) { // obtain the sums required for the posterior parameters
    
    uword patt_i = id_patt.at(i); // id missing outcome pattern
    
    if (i < patt_count) { // obtain all unique vcov_inv matrices required for the sums in the posterior parameters
      
      mat L_patt_inv = inv( trimatl( chol_update(L_D, ind_RE_patt.at(i)) ) );
      D_inv.at(i) =  L_patt_inv.t() * L_patt_inv;
      
    }

    mat X_dot_i = X_dot.submat(ind_RE_patt.at(patt_i) - 1 + i*re_count, 
                               ind_FE_patt.at(patt_i) - 1); 
    
    vec u_i = b_mat.row(i).t() + X_dot_i * current_FE_in_hc.elem(ind_FE_patt.at(patt_i) - 1);
    mat J_i = J.cols(ind_FE_patt.at(patt_i) - 1);
    D_inv_i = D_inv.at(patt_i);  
    
    mat XD_i = X_dot_i.t() * D_inv_i;
    mat XDX_i = XD_i * X_dot_i;
    
    vec sum_JXDu  += add_zero_colrows(XD_i, p_hc, XD_i.n_cols, ind_FE_patt.at(patt_i) - 1, regspace<uvec>(0,  XD_i.n_cols - 1)) * u_i; // add zero-rows
    mat sum_JXDXJ += add_zero_colrows(XDX_i, p_hc, p_hc, ind_FE_patt.at(patt_i) - 1, ind_FE_patt.at(patt_i) - 1); // add zero-rows and zero-columns
    
  }
  
  vec mean_1; // posterior mean vector
  mat Tau_1; // posterior vcov matrix
  
  Tau_1 = (prior_Tau_fe_in_hc + sum_JXDXJ).i();
  mean_1 = Tau_1 * (prior_Tau_fe_in_hc * prior_mean_fe_in_hc + sum_JXDu);
  cube U_1(p_hc, p_hc, 1); U_1.slice(0) = chol(Tau_1); // propose_mvnorm_mat() expects a cube
  
  res_betas.submat(it, ind_FE_in_hc - 1) = (propose_mvnorm_mat(1, U_1, {1}) + mean_1).t();
  
  // FE not in HC - Metropolis-Hastings sampling
  
  res_betas.submat(it, ind_FE_notin_hc - 1) = 99; //! <---------------------------------
  
  //
  
  betas = mat2field_mat(res_betas.row(it), ind_FE - 1); //! check if this function works with field of vectors <---------------------------------
  
)

#endif

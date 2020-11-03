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
                   const vec &prior_mean_FE, const mat &prior_invSigma_FE, // prior
                   const mat &b_mat; // it-th sampled RE
                   const mat &L, // RE corr matrix factorization
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
  uword n = b_mat.n_rows(); // number of unique subjects
  uword re_count = L.n_cols(); // number of random effects
  uword patt_count = id_patt.max(); // number of unique outcome-missing patterns          
  vec prior_mean_FE_in_hc = prior_mean_FE.elem(ind_FE_in_hc - 1); // prior fot FE in HC 
  mat prior_invSigma_FE_in_hc = prior_invSigma_FE.submat(ind_FE_in_hc - 1, ind_FE_in_hc - 1);
  vec prior_mean_FE_notin_hc = prior_mean_FE.elem(ind_FE_notin_hc - 1); // prior for FE NOT in HC 
  mat prior_invSigma_FE_notin_hc = prior_invSigma_FE.submat(ind_FE_notin_hc - 1, ind_FE_notin_hc - 1);
  
  // FE in HC - Gibss sampling
  uword p_hc = ind_FE_in_hc.n_elem(); // number of FE in HC
  mat::fixed<p_hc, p_hc> J; J.eye(); // identity matrix, required in JXDXJ_i and JXDu_i
  mat::fixed<p_hc, p_hc> sum_JXDXJ; sum_JXDXJ.zeros(); // sum required for posterior parameters
  vec::fixed<p_hc> sum_JXDu; sum_JXDu.zeros(); // sum required for posterior parameters
  
  
  mat::fixed<re_count, re_count> L_D = L.each_row() % sds.t(); // RE vcov matrix factorization
  mat::fixed<re_count, re_count> D = L_D * L_D.t(); // RE vcov matrix
  field<mat> D_inv(patt_count, 1); // all unique vcov_inv matrices accross the missing outcome patterns

  for (uword i = 0; i < n; ++i) { // obtain the sums required for the posterior parameters
    
    uword patt_i = id_patt.at(i); // id missing outcome pattern
    
    if (i < patt_count) { // obtain all unique vcov_inv matrices required for the sums in the posterior parameters
      
      D_inv(i, 0) = D.submat( ind_RE_patt(i, 0) - 1, ind_RE_patt(i, 0) - 1).i();
      
    }

    mat X_dot_i = X_dot.submat(ind_RE_patt(patt_i, 0) - 1 + i*re_count, 
                               ind_FE_patt(patt_i, 0) - 1); 
    vec::fixed<re_count> u_i = b_mat.row(i); //! b is not u <---------------------------------
    mat J_i = J.cols(ind_FE_patt(patt_i, 0) - 1);
    D_inv_i = D_inv(patt_i, 0);  
    
    mat sum_JXDu += J_i * (X_dot_i.t() * D_inv_i * X_dot_i) * J_i.t();
    sum_JXDXJ    += J_i * (X_dot_i.t() * D_inv_i) * u_i;
    
  }
  
  mat::fixed<1, p_hc> mean_1; // posterior mean vector
  mat::fixed<p_hc, p_hc, 1> Tau_1; // posterior vcov matrix (propose_mvnorm_mat() expects a cube)
  
  Tau_1 = (prior_invTau_fe_in_hc + sum_JXDXJ).i();
  mean_1 = Tau_1 * (prior_invTau_fe_in_hc * prior_mean_fe_in_hc + sum_JXDu);
  L_1 = chol(Tau_1);
  
  res_betas.submat(it, ind_FE_in_hc - 1) = propose_mvnorm_mat(1, L_1, {1}).t();
  
  // FE not in HC - Metropolis-Hastings sampling
  
  res_betas.submat(it, ind_FE_notin_hc - 1) = 99; //! <---------------------------------
  
  //
  
  betas = mat2field_mat(res_betas.row(it), ind_FE - 1); //! check if this function works with field of vectors <---------------------------------
  
)

#endif

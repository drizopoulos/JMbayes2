#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

field<uvec> List2Field_uvec (const List &uVecs, bool substract1 = true) {
  uword n_list = uVecs.size();
  field<uvec> res(n_list);
  for (uword i = 0; i < n_list; ++i) {
    if (substract1) {
      res.at(i) = as<arma::uvec>(uVecs[i]) - 1;
    } else {
      res.at(i) = as<arma::uvec>(uVecs[i]);
    }
  }
  return res;
}

mat cov2cor (const mat &V) {
  vec Is = sqrt(1.0 / V.diag());
  mat out = V.each_col() % Is;
  out.each_row() %= Is.t();
  return out;
}

mat docall_cbindF (const field<mat> &Mats) {
  uword n = Mats.n_elem;
  uvec ncols(n);
  for (uword k = 0; k < n; ++k) {
    ncols.at(k) = Mats.at(k).n_cols;
  }
  uword N = sum(ncols);
  uword col_start = 0;
  uword col_end = ncols.at(0) - 1;
  mat out(Mats.at(0).n_rows, N);
  for (uword k = 0; k < n; ++k) {
    if (k > 0) {
      col_start += ncols.at(k - 1);
      col_end += ncols.at(k);
    }
    out.cols(col_start, col_end) = Mats.at(k);
  }
  return out;
}

field<mat> List2Field_mat (const List &Mats) {
  uword n_list = Mats.size();
  field<mat> res(n_list);
  for (uword i = 0; i < n_list; ++i) {
    res.at(i) = as<mat>(Mats[i]);
  }
  return res;
}

field<vec> List2Field_vec (const List &Vecs) {
  uword n_list = Vecs.size();
  field<vec> res(n_list);
  for (uword i = 0; i < n_list; ++i) {
    res.at(i) = as<vec>(Vecs[i]);
  }
  return res;
}

field<vec> vec2field (const vec &betas, const field<uvec> &ind_FE) {
  uword n = ind_FE.n_elem;
  field<vec> out(n);
  for (uword i = 0; i < n; i++) {
    out.at(i) = betas.rows(ind_FE.at(i));
  }
  return out;
}

mat rnorm_mat (const uword& rows, const uword& cols) {
  mat out(rows, cols);
  out.each_col([&](vec& x) {x = as<vec>(rnorm(rows)); } );
  return out;
}

mat propose_mvnorm_vec (const int &n, // number of samples
                        const mat &U, // upper triangular Choleski factorization of the vcov matrix
                        const double &scale) {
  uword ncols = U.n_cols;
  mat res(n, ncols);
  res = scale * (rnorm_mat(n, ncols) * U);
  return res.t(); // j-th column reports the j-th sample
}

vec docall_rbindF (const field<vec> &F) { // binds a field of vectors into one vector
  uword n = F.n_elem;
  uword nrows = 0;
  uvec rows(n);
  for (uword i = 0; i < n; i++) {
    rows.at(i) = F.at(i).n_rows;
    nrows += rows.at(i);
  }
  vec V(nrows);
  uword ii = 0;
  for (uword i = 0; i < n; i++) {
    V.rows(ii, ii - 1 + rows.at(i)) = F.at(i);
    ii += rows.at(i);
  }
  return V;
}

// ?? The function below could be optimizaed to update in blocks, rather then element by element
mat add_zero_colrows (const mat &M, // adds zero-rows and/or zero-cols to a matrix M
                      const uword &nrows, // n_rows in the target matrix
                      const uword &ncols, // n_cols in the target matrix
                      const uvec &rows_ind,    // ind where to place the M's rows (zero-rows will be 'added' to the absent ind). the number of ind must match the M's n_rows
                      const uvec &cols_ind) { // ind where to place the M's cols (zero-cols will be 'added' to the absent ind). the number of ind must match the M's n_cols
  
  mat Res(nrows, ncols, fill::zeros);
  uword M_nrows = M.n_rows;
  uword M_ncols = M.n_cols;
  
  for(uword i = 0; i < M_nrows; i++) { // by row
    for(uword j = 0; j < M_ncols; j++) { // by col
      Res.at(rows_ind.at(i), cols_ind.at(j)) = M.at(i, j);
    }
  }
  return Res;
}

mat add_zero_rows (const mat &M, // adds zero-rows to a matrix M
                   const uword &nrows, // n_rows in the target matrix
                   const uvec &rows_ind) { // ind where to place the M's rows (zero-rows will be 'added' to the absent ind). the number of ind must match the M's n_rows
  
  mat Res(nrows, M.n_cols, fill::zeros);
  uword M_nrows = M.n_rows;
  
  for(uword j = 0; j < M_nrows; j++) { // by col
    Res.row(rows_ind.at(j)) = M.row(j);
  }
  return Res;
}

mat rank1_update(const mat &M, // performs rank-1 update: Res = M + v * u.t() 
                 const vec &v) {
  
  uword n = v.n_elem;
  mat Res = M;
  vec v2   = v;
  
  for (uword i = 0; i < n; i++) {
    double r = pow( pow(Res.at(i, i), 2) + pow(v2.at(i), 2), 0.5);
    double c = r / Res.at(i, i);
    double s = v2.at(i) / Res.at(i, i);
    Res.at(i, i) = r;
    
    if (i < n - 1) {
      Res.submat(i + 1, i, n - 1, i) = (Res.submat(i + 1, i, n - 1, i) + s * v2.rows(i + 1, n - 1)) / c;
      v2.rows(i + 1, n - 1) = c * v2.rows(i + 1, n - 1) - s * Res.submat(i + 1, i, n - 1, i);
    }
  }
  return Res;
}

mat chol_update(const mat &L, // If L = chol(M), returns chol(M.submat(keep, keep))
                const uvec &keep) { // keep must be a sorted vector, i.e, {2, 4, 5}, and counts from 0
  
  // later we can try to extend this approach further to obtain inv(L) from the required inv(L_i)
  
  uvec rem = regspace<uvec>(0,  L.n_cols - 1); rem.shed_rows(keep); // cols-rows to remove
  mat Res = L;
  uword n = rem.n_elem;
  
  for (uword i = 0; i < n; i++) { // rank-1 update for each col-row to be removed
    
    uword last_col = Res.n_cols - 1;
    
    if(rem.at(i) < last_col) {
      Res.submat(rem.at(i) + 1, rem.at(i) + 1, last_col, last_col) = rank1_update(Res.submat(rem.at(i) + 1, rem.at(i) + 1, last_col, last_col),
                 Res.submat(rem.at(i) + 1, rem.at(i), last_col, rem.at(i)));
    }
    
    Res.shed_row(rem.at(i));
    Res.shed_col(rem.at(i));
    rem = rem - 1;
  }
  return Res;
}


////////////////////////////////////////////////////////////////////////////////


void update_betas (field<vec> &betas, // it-th sampled fixed effects
                   mat &res_betas, // all sampled fixed effects
                   const uword &it, // current iteration
                   const vec &prior_mean_betas_HC, 
                   const mat &prior_Tau_betas_HC, 
                   const mat &b_mat, // it-th sampled RE
                   const mat &L, // RE corr matrix factorization (upper matrix)
                   const vec &sds, // RE SDs
                   const mat &X_dot, // X_dot matrix
                   const field<uvec> &ind_FE, // indices for the FE in res_betas[it,] belonging to the field betas. E.g., {{1,2,3}, {4, 5}, {6}}
                   const uvec &ind_FE_HC, // indices for the FE present in the HC (cols in res_betas). Counts from 0.
                   const uvec &id_patt, // vector with the ids' outcome missing pattern
                   const field<uvec> &ind_RE_patt, // indices for the RE present in each outcome missing pattern (cols in D). Counts from 0.
                   const field<uvec> &ind_FE_patt // indices for the FE (in HC) present in each outcome missing pattern (cols in X_dot). Counts from 0.
) {
  
  vec betas_vec = docall_rbindF(betas);
  
  // FE in HC - Gibss sampling
  
  uword n = b_mat.n_rows; // number of unique subjects
  uword patt_count = id_patt.max() + 1; // number of unique outcome-missing patterns          
  uword p_HC = ind_FE_HC.n_elem; // number of HC-FE
  uword q = b_mat.n_cols; // number of random effects
  mat sum_JXDXJ(p_HC, p_HC, fill::zeros); // sum for the posterior parameters
  vec sum_JXDu(p_HC, fill::zeros); // sum for the posterior parameters
  
  mat U = L.each_row() % sds.t(); // RE vcov matrix Choesky factorization (upper) //?? check this
  
  field<mat> D_inv(patt_count); // all unique vcov_inv matrices accross the missing outcome patterns

  for (uword i = 0; i < n; ++i) { // i-th patient
    
    uword patt_i = id_patt.at(i); // id missing outcome pattern // uword
    uvec ind_FE_i = ind_FE_patt.at(patt_i); // uvec
    uvec ind_RE_i = ind_RE_patt.at(patt_i);
    
    if (i < patt_count) { // obtain all unique vcov_inv matrices required for the sums in the posterior parameters
      
      mat U_patt_inv = inv( trimatu( chol_update(U, ind_RE_patt.at(i)) ) ); // mat
      D_inv.at(i) =  U_patt_inv * U_patt_inv.t(); // mat //?? do you know a better way to do this caculation?
      
    }
    
    mat X_dot_i = X_dot.rows(i*q, (i + 1)*q - 1); // mat
    X_dot_i = X_dot_i.submat(ind_RE_i, ind_FE_i);
    vec b_i = b_mat.row(i).t();
    vec u_i = b_i.rows(ind_RE_i) + X_dot_i * betas_vec.rows(ind_FE_i);
    mat D_inv_i = D_inv.at(patt_i); // mat 
    mat XD_i = X_dot_i.t() * D_inv_i; // mat
    mat XDX_i = XD_i * X_dot_i; // mat
    sum_JXDu  += add_zero_rows(XD_i*u_i, p_HC, ind_FE_i); // mat
    sum_JXDXJ += add_zero_colrows(XDX_i, p_HC, p_HC, ind_FE_i, ind_FE_i); // mat
  }
  
  mat Tau_1 = inv(prior_Tau_betas_HC + sum_JXDXJ); //?? Can I write the Tau_1 in terms of an L matrix? to avoid the use of chol(Tau_1) later <---------------------------------
  vec mean_1 = Tau_1 * (prior_Tau_betas_HC * prior_mean_betas_HC + sum_JXDu);
  mat U_1 = chol(Tau_1);
  betas_vec.rows(ind_FE_HC) = propose_mvnorm_vec(1, U_1, 1) + mean_1; // vec
  betas = vec2field(betas_vec, ind_FE); // field<vec>  
  res_betas.row(it) = betas_vec.t(); // vec.t()
  
}

// [[Rcpp::export]]

List gibbs_sim (List data) {
  

  field<vec> betas = List2Field_vec(as<List>(data["betas"]));
  vec betas_vec = docall_rbindF(betas);
  uword n_betas = betas_vec.n_rows;
  uword n_iter = as<uword>(data["n_iter"]);
  mat res_betas(n_iter, n_betas, fill::zeros);
  vec prior_mean_betas_HC = as<vec>(data["mean_betas_HC"]);
  mat prior_Tau_betas_HC = as<mat>(data["Tau_betas_HC"]);
  field<mat> b = List2Field_mat(as<List>(data["b"]));
  mat b_mat = docall_cbindF(b);
  mat D = as<mat>(data["D"]);
  vec sds = sqrt(D.diag());
  mat R = cov2cor(D);
  mat L = chol(R);
  mat X_dot = as<mat>(data["X_dot"]);
  field<uvec> ind_FE = List2Field_uvec(as<List>(data["ind_FE"]), true);
  uvec ind_FE_HC = as<uvec>(data["ind_FE_HC"]) - 1;
  uvec id_patt = as<uvec>(data["id_patt"]) - 1;
  field<uvec> ind_RE_patt = List2Field_uvec(as<List>(data["ind_RE_patt"]), true);
  field<uvec> ind_FE_patt = List2Field_uvec(as<List>(data["ind_FE_patt"]), true);
  
  for (uword i = 0; i < n_iter; ++i) {
    
    update_betas(betas,
                 res_betas,
                 i,
                 prior_mean_betas_HC, 
                 prior_Tau_betas_HC, 
                 b_mat, 
                 L,
                 sds,
                 X_dot,
                 ind_FE,
                 ind_FE_HC,
                 id_patt,
                 ind_RE_patt,
                 ind_FE_patt);
  }
  
  return List::create(
    Named("betas") = betas,
    Named("res_betas") = res_betas
  );
}


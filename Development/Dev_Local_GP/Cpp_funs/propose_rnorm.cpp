//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

mat rnorm_mat (const uword& rows, const uword& cols) {
  mat out(rows, cols);
  out.each_col([&](vec& x) {x = as<vec>(rnorm(rows)); } );
  return out;
}

// [[Rcpp::export]]
mat propose_mvnorm_cube (const uword& n, const cube& S, const vec& sigmas) {
  int ncol_per_slice = S.n_cols;
  int slices = S.n_slices;
  cube out(n, ncol_per_slice, slices);
  for (uword i = 0; i < slices; i++) {
    out.slice(i) = sqrt(sigmas.at(i)) * (rnorm_mat(n, ncol_per_slice) * chol(S.slice(i)));
  }
  mat outnew = out.row(0);
  //cube Z(1, ncol_per_slice, slices, fill::randn);
  //cube out = S; 
  //out.each_slice([&] (mat& X) {mat Z(1, ncol_per_slice, fill::randn); 
  //  X = Z * trans(chol(X)); } );
  
  
  //cube S_transchol = S;
  //S_transchol.each_slice([] (mat& X) {X = trans(chol(X)); } );
  //.each_slice( [](mat& A) {A = A + 1; } );
  //Rcout << S.slice(0) << endl;
  //Rcout << S_transchol.slice(0) << endl;
  return outnew.t();
}
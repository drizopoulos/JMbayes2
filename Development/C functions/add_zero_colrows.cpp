#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// Tried to develop a new version for the add_zero_colrows() with block updating, but it seems to not be better than using nested for-cycles

// [[Rcpp::export]]
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

// [[Rcpp::export]]
mat add_zero_colrows2 (const mat &M, // adds zero-rows and/or zero-cols to a matrix M
                      const uword &nrows, // n_rows in the target matrix
                      const uword &ncols, // n_cols in the target matrix
                      const uvec &rows_ind,    // ind where to place the M's rows (zero-rows will be 'added' to the absent ind). the number of ind must match the M's n_rows
                      const uvec &cols_ind) { // ind where to place the M's cols (zero-cols will be 'added' to the absent ind). the number of ind must match the M's n_cols
  mat Res(nrows, ncols, fill::zeros);
  Res(rows_ind, cols_ind) = M;
  return Res;
}

// [[Rcpp::export]]
mat add_zero_colrows3 (const mat &M, // adds zero-rows and/or zero-cols to a matrix M
                       const uword &nrows, // n_rows in the target matrix
                       const uword &ncols, // n_cols in the target matrix
                       const uvec &rows_ind,    // ind where to place the M's rows (zero-rows will be 'added' to the absent ind). the number of ind must match the M's n_rows
                       const uvec &cols_ind) { // ind where to place the M's cols (zero-cols will be 'added' to the absent ind). the number of ind must match the M's n_cols
  mat Res(nrows, ncols, fill::zeros);
  Res.submat(rows_ind, cols_ind) = M;
  return Res;
}

/*** R
M <- matrix(runif(10000), 100, 100) # initial matrix 
nrows <- 200 # n_rows in the target matrix
ncols <- 200 # n_cols in the target matrix
rows_ind <- sample(1:200, 100) # rows ind where to place the M's rows (zero-rows will be 'added' to the absent ind). the number of ind must match the M's n_rows
cols_ind <- sample(1:200, 100) # cols ind where to place the M's rows (zero-rows will be 'added' to the absent ind). the number of ind must match the M's n_cols

M1 <- add_zero_colrows(M, nrows, ncols, rows_ind - 1, cols_ind - 1)
M2 <- add_zero_colrows2(M, nrows, ncols, rows_ind - 1, cols_ind - 1)
M3 <- add_zero_colrows3(M, nrows, ncols, rows_ind - 1, cols_ind - 1)
all(M1 == M2)
all(M1 == M3)

set.seed(2022)
microbenchmark::microbenchmark(
  add_zero_colrows(M, nrows, ncols, rows_ind - 1, cols_ind - 1),
  add_zero_colrows2(M, nrows, ncols, rows_ind - 1, cols_ind - 1),
  add_zero_colrows3(M, nrows, ncols, rows_ind - 1, cols_ind - 1),
  times = 100000L
)
*/

#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

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
mat add_zero_rows (const mat &M, // adds zero-rows to a matrix M
                   const uword &nrows, // n_rows in the target matrix
                   const uvec &rows_ind) { // ind where to place the M's rows (zero-rows will be 'added' to the absent ind). the number of ind must match the M's n_rows
  
  mat Res(nrows, M.n_cols, fill::zeros);
  uword M_nrows = M.n_rows;
  
  for(uword j = 0; j < M_nrows; j++) {
    Res.row(rows_ind.at(j)) = M.row(j);
  }
  return Res;
}


/*** R
nrows <- 3

M <- matrix(99, nrow = nrows, ncol = nrows)

all_rows_ind <- lapply(seq_len(nrows), combn, x = seq_len(nrows), simplify = FALSE)
all_rows_ind <- unlist(all_rows_ind, recursive = FALSE)


# add_zero_rows
for(rows_ind in all_rows_ind) {
  
  if(length(rows_ind) == nrows) next
  
  R <- M
  R[rows_ind, ] <- 0
  M_sub <- M[-rows_ind, , drop = FALSE]
  
  C <- add_zero_rows(M_sub, nrows, setdiff(seq_len(nrows), rows_ind)-1)
  
  cat(paste0(all.equal(C, R), " "))
  if(!all.equal(C, R)) print(paste0(rows_ind, "\n"))
  
}


# add_zero_colrows
for(rows_ind in all_rows_ind) {
  
  if(length(rows_ind) == nrows) next
  
  R <- M
  R[rows_ind,] <- 0
  R[,rows_ind] <- 0
  M_sub <- M[-rows_ind, -rows_ind, drop = FALSE]
  
  C <- add_zero_colrows(M_sub, nrows, nrows, setdiff(seq_len(nrows), rows_ind)-1, setdiff(seq_len(nrows), rows_ind)-1)

  cat(paste0(all.equal(C, R), " "))
  if(!all.equal(C, R)) print(paste0(rows_ind, "\n"))
  
}

*/

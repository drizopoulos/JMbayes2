#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List test_lists(List L, unsigned int elem) {
    List L_i = L[elem];
    return L_i;
}



/*** R
l <- list("a", list("b", "c"))

test_lists(l, 0)

test_lists(l, 1)
*/

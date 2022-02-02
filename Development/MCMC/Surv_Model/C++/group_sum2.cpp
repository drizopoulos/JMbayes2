#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
vec group_sum2 (const vec &x, const uvec &ind) {
    uvec tt = unique(ind);
    uword n = tt.n_rows;
    uvec ind2 = ind - 1;
    vec out(n);
    for (uword i = 0; i < n; ++i) {
        out.at(i) = sum(x.elem(find(ind2 == i)));
    }
    return out;
}



/*** R
xx <- rnorm(30)
ind <- sample(1:5, 30, TRUE)

rowsum(xx, ind)
group_sum2(xx, ind)
*/

library(Rcpp)
library(RcppArmadillo)

sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\arma_logmat.cpp')

x <- matrix(1:4, ncol = 2, nrow = 2)

arma_logmat(x, 0, 0)

arma_logmat(x, 1, 1)
log(1)
log(2)
log(3)
log(4)

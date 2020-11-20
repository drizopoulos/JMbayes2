library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(ggplot2)

sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\docall_cbind.cpp')

listmats <- list(matrix(rnorm(200), nrow = 100, ncol = 2), 
                 matrix(rnorm(100, 100, 0.5), nrow = 100, ncol = 1))

docall_cbindL(listmats)

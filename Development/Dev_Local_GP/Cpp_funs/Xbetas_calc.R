library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(ggplot2)

sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\Xbetas_calc.cpp')

klain2 <- Xbeta_calc(X, betas)

all.equal(klain1, klain2[[1]])

klain1 <- X[[1]] %*% betas[[1]]
class(klain1)

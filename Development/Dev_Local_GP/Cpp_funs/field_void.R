library(Rcpp)
library(RcppArmadillo)

sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\field_void.cpp')

X <- list(1:9, 100:200)

res <- test_void(X)

res[1, 1]
res[2, 1]

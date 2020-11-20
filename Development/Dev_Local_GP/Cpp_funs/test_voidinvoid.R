library(Rcpp)
library(RcppArmadillo)

sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/voidinvoid.cpp')

X <- 1:4

test_void(X)

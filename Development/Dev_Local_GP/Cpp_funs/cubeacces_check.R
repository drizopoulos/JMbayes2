library(Rcpp)
library(RcppArmadillo)
sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/cubeacces.cpp')

res <- matrix(1:40, nrow = 10, ncol = 4)

chk <- foo(res)

chk[, ,1]

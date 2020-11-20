JMbayes2:::docall_cbind(obj$initial_values$b)

klain2 <- array(0.0, dim = c(5000, 493, 4))

dim(klain2[, , ])
dim(res)
dim(klain2)

library(Rcpp)
library(RcppArmadillo)
sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/cubeacces.cpp')


res <- matrix(1:40, nrow = 10, ncol = 4)
res[4] <- 10

foo(klain2, res)
klain3 <- foo(klain2, res)
dim(klain2)
dim(klain3)

klain3[1, ,1]
?array

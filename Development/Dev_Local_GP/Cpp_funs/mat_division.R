library(Rcpp)
library(RcppArmadillo)

x <- c(3, 4, 6, 7, 1, 5)

sd(x)

mean(sum(x))
mean(x)

xmat <- matrix(10.0, nrow = 200, ncol = 4)
xmat / 10

sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/mat_division.cpp')


foo(xmat)

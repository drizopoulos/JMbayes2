library(Rcpp)
library(RcppArmadillo)

sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\test_if_for_scope.cpp')

x <- list(1:4, 5:8, 9:12)
chk <- foo(x, 1)
chk[[1]]

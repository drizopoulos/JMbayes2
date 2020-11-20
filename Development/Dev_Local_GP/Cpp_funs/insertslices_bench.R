library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(ggplot2)

sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/insertslices.cpp')
sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/insertslices_2.cpp')
sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/insertslices_3.cpp')

x <- array(0.0, dim = c(100, 4, 1))
x2 <- foo(x)

dim(x)
dim(x2)


bench <- microbenchmark('foo' = {foo()}, 
                        'foo2' = {foo2()}, 
                        'foo3' = {foo3()})

bench
autoplot(bench)

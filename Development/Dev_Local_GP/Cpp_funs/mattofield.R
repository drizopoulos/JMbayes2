library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(ggplot2)

sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\mat_to_field.cpp')

indRE <- list(1:2, 3:4, 5, 6)
indRE <- lapply(indRE, FUN = function (x) x - 1)

klain4 <- mat_to_field(chk1, indRE)

klain4[[1]]

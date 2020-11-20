library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(ggplot2)

sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\rnorm_R.cpp')
sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\rnorm2_R.cpp')
sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\rnorm3_R.cpp')


rnorm(5, 5, 0, 1)
rnorm2(5, 5, 0, 1)
rnorm3(5, 5)


system.time({
  rnorm(5000000, 1000, 0, 1)
})

system.time({
  rnorm2(5000000, 1000, 0, 1)
})

bench <- microbenchmark('rnorm' = {rnorm(100000, 100, 0, 1)}, 
          'rnorm2' = {rnorm2(100000, 100, 0, 1)}, 
          'rnorm3' = {rnorm3(100000, 100)}
)

bench
autoplot(bench)
save(bench, 
     file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\bench_rnorm.RData')

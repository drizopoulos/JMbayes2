library(Rcpp)
library(RcppArmadillo)

sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/iter_mean/test_basic.cpp')
sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/iter_mean/test_iter.cpp')

testobj <- list(matrix(rnorm(300), nrow = 100, ncol = 3), matrix(rnorm(300), nrow = 100, ncol = 3), 
     matrix(rnorm(300), nrow = 100, ncol = 3), matrix(rnorm(300), nrow = 100, ncol = 3))

chk1 <- foo(testobj)
chk2 <- foo(testobj)
chk3 <- (testobj[[1]] + testobj[[2]] + testobj[[3]] + testobj[[4]]) / 4

all.equal(chk1, chk2)
all.equal(chk1, chk3)
all.equal(chk2, chk3)


out = testobj[[1]]
for (i in 2:4) {
  out = (testobj[[i]] - out)
}

all.equal(chk2[1, ], out[1, ])


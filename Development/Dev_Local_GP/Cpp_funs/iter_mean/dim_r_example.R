library(Rcpp)
library(RcppArmadillo)

sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_Funs/iter_mean/foo_outer.cpp')


x <- matrix(rnorm(5 * 100), 100, 5)
K <- 9
cum_sum <- numeric(5)
out_prod <- matrix(0.0, 5, 5)
for (i in 1:K) {
  cum_sum <- cum_sum + x[i, ]
  out_prod <- out_prod + x[i, ] %o% x[i, ]
}
means <- cum_sum / K
dim((out_prod / K - means %o% means) * K / (K - 1))
var(x[1:K, ])

x[9, ]

chk <- x[1, ] %o% x[1, ]
chk3 <- kronecker(x[1, ], x[1, ])
chk2 <- foo(x[1, ])

all.equal(chk, chk2)


x[1, ] %o% x[1, ]
foo(x[1, ])


x <- list(matrix(rnorm(5 * 100), 100, 5), matrix(rnorm(5 * 100), 100, 5), matrix(rnorm(5 * 100), 100, 5))


cum_sum <- x[[1]]
out_prod <- x[[1]] %o% x[[1]]

for(i in 2:3) {
  cum_sum = cum_sum + x[[i]]
  out_prod = out_prod + x[[i]] %o% x[[i]]
}

out_prod[1, , 1, ]

klain1 <- x[[1]] %o% x[[1]]
klain1[1, ,1 , ]
x[[1]][1, ] %o% x[[1]][1, ]
class(klain1)
klain2 <- foo(x[[1]])
class(klain2)
dim(klain2)

means = cum_sum / 3
dim((out_prod / 3 - means %o% means) * K / (K - 1))

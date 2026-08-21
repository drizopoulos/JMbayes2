k <- 6; n <- 2500
xx <- matrix(rnorm(n * k), n, k)
VV <- var(xx)
LL <- chol(VV)
yy <- rnorm(k)
kk <- c(3, 5)


tt1 <- chol_update(LL, kk)
tt2 <- chol_update_fast(LL, kk)

all.equal(tt1, tt2)

library("rbenchmark")

benchmark(
    Cpp1 = chol_update(LL, kk),
    Cpp2 = chol_update_fast(LL, kk),
    replications = 20000
)






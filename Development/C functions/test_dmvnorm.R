k <- 5; n <- 2500
xx <- matrix(rnorm(n * k), n, k)
VV <- var(xx)
LL <- chol(VV)

n <- 1200
k <- 15
bb <- matrix(rnorm(n * k), n, k)
sds <- abs(rnorm(k))

tt1 <- scaled_Armadillo(bb, sds)
tt2 <- scaled_pointers(bb, sds)

all.equal(tt1, tt2)

library("rbenchmark")

benchmark(
    Armadillo = scaled_Armadillo(bb, sds),
    Pointers = scaled_pointers(bb, sds),
    replications = 50000
)






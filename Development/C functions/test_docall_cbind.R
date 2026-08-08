k <- 20; n <- 2900
xx <- lapply(1:10, function (i) matrix(rnorm(n * k), n, k))


tt1 <- docall_cbindF(xx)
tt2 <- docall_cbindF_fast(xx)

all.equal(tt1, tt2)

library("rbenchmark")

benchmark(
    Cpp1 = docall_cbindF(xx),
    Cpp2 = docall_cbindF_fast(xx),
    replications = 1000
)






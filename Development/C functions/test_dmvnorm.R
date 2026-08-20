k <- 5; n <- 2500
xx <- matrix(rnorm(n * k), n, k)
VV <- var(xx)
LL <- chol(VV)


tt1 <- log_dmvnrm_chol(xx, LL)
tt2 <- log_dmvnrm_chol_optimal(xx, LL)

all.equal(tt1, tt2)

library("rbenchmark")

benchmark(
    Cpp1 = log_dmvnrm_chol(xx, LL),
    Cpp2 = log_dmvnrm_chol_optimal(xx, LL),
    replications = 20000
)






xx <- rnorm(1e6)
mm <- rnorm(1e6)
sig <- 0.5

tt1 <- log_dnorm(xx, mm, sig)
tt2 <- log_dnorm_fast(xx, mm, sig)

all.equal(tt1, tt2)

library("rbenchmark")

benchmark(
    Cpp1 = log_dnorm(xx, mm, sig),
    Cpp2 = log_dnorm_fast(xx, mm, sig),
    replications = 1000
)

xx <- rnorm(1e5)
lnk <- "identity"

tt1 <- mu_fun(xx, lnk)
tt2 <- mu_fun_fast(xx, lnk)

all.equal(tt1, tt2)

library("rbenchmark")

benchmark(
    Cpp1 = mu_fun(xx, lnk),
    Cpp2 = mu_fun_fast(xx, lnk),
    replications = 1000
)






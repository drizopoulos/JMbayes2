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

#########################

nn <- 1e5
prs <- runif(nn)
szs <- sample(1:3000, nn, TRUE)
xx <- rbinom(nn, szs, prs)

tt1 <- log_dbinom(xx, szs, prs)
tt2 <- log_dbinom_ultra(xx, szs, prs)

all.equal(tt1, tt2)

library("rbenchmark")

benchmark(
    Cpp1 = log_dbinom(xx, szs, prs),
    Cpp2 = log_dbinom_ultra(xx, szs, prs),
    replications = 1000
)

#########################

nn <- 1e5
prs <- runif(nn)
phi <- 2
szs <- sample(1:3000, nn, TRUE)
xx <- rbinom(nn, szs, prs)

tt1 <- log_dbbinom(xx, szs, prs, phi)
tt2 <- log_dbbinom_ultra(xx, szs, prs, phi)

all.equal(tt1, tt2)

library("rbenchmark")

benchmark(
    Cpp1 = log_dbbinom(xx, szs, prs, phi),
    Cpp2 = log_dbbinom_ultra(xx, szs, prs, phi),
    replications = 100
)






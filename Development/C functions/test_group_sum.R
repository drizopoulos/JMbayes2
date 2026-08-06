xx <- rnorm(1e6)
ii <- sort(sample(1:500, 1e6, TRUE))
jj <- create_fast_ind(ii - 1)

v1 <- c(rowsum(xx, ii, reorder = FALSE))
v2 <- c(group_sum(xx, jj))
v3 <- c(group_sum2(xx, ii))
v4 <- c(group_sum2(xx, jj))

all.equal(v1, v2)
all.equal(v1, v3)
all.equal(v1, v4)

library("rbenchmark")

benchmark(
    R = c(rowsum(xx, ii, reorder = FALSE)),
    Cpp1 = c(group_sum(xx, jj)),
    Cpp2 = c(group_sum2(xx, jj)),
    replications = 1000
)

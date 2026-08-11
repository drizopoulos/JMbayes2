v1 <- robbins_monro(0.9, 1, 20)
v2 <- robbins_monro_fast(0.9, 1, 20)

all.equal(v1, v2)

library("rbenchmark")

benchmark(
    Cpp1 = robbins_monro(0.5, 1, 2),
    Cpp2 = robbins_monro_fast(0.5, 1, 2),
    replications = 1000
)

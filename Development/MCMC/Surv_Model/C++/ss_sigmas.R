log_sigmas <- c(0.5, -20, -0.8, -20.0, 0)
has_sigmas <- as.integer(log_sigmas > -20)
ss_sigmas <- c(TRUE, FALSE, FALSE, FALSE, TRUE)
idL <- replicate(5, rep(1:20, sample(1:20, 20)), simplify = FALSE)
ns <- lapply(idL, function (x) length(unique(x)))
Ns <- lapply(idL, length)
#########################
sigmas <- split(exp(log_sigmas), seq_along(log_sigmas))
sigmas[ss_sigmas] <- mapply(rep, x = sigmas[ss_sigmas],
                            length.out = ns[ss_sigmas], SIMPLIFY = FALSE)
sigmas[] <- lapply(sigmas, jitter)

create_sigma_list <- function (sigmas, ss_sigmas, idL) {
    n <- length(sigmas)
    out <- vector("list", n)
    for (i in 1:n) {
        sigmas_i <- sigmas[[i]]
        id_i <- idL[[i]]
        if (ss_sigmas[i]) {
            out[[i]] <- sigmas_i[id_i]
        } else {
            out[[i]] <- rep(sigmas_i, length(id_i))
        }
    }
    out
}

SS1 <- create_sigma_list(sigmas, ss_sigmas, idL)
SS2 <- create_sigmas_field(sigmas, ss_sigmas, lapply(idL, function (x) x - 1))

all.equal(c(SS2[1, ][[1]]), SS1[[1]])
all.equal(c(SS2[2, ][[1]]), SS1[[2]])
all.equal(c(SS2[3, ][[1]]), SS1[[3]])
all.equal(c(SS2[4, ][[1]]), SS1[[4]])
all.equal(c(SS2[5, ][[1]]), SS1[[5]])


library("rbenchmark")
idL2 <- lapply(idL, function (x) x - 1)
benchmark(R = create_sigma_list(sigmas, ss_sigmas, idL),
          Cpp = create_sigmas_field(sigmas, ss_sigmas, idL2),
          replications = 120000)

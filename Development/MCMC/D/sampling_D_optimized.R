library("survival")
library("nlme")
library("GLMMadaptive")
library("splines")
data("pbc2", package = "JM")
data("pbc2.id", package = "JM")
source(file.path(getwd(), "R/jm.R"))
source(file.path(getwd(), "R/help_functions.R"))
source(file.path(getwd(), "Development/jm/R_to_Cpp.R"))
source(file.path(getwd(), "Development/jm/PBC_data.R"))
source(file.path(getwd(), "Development/MCMC/D/sampling_D_Funs.R"))

source(file.path(getwd(), "Development/MCMC/D/D_examples.R"))

##########################################################################################

# Select D matrix from the D_examples
D <- D4

p <- ncol(D)
K <- as.integer(round(p * (p - 1) / 2))
sds <- sqrt(diag(D))
R <- cov2cor(D)
inv_R <- solve(R)
L <- chol(R)
Ip <- diag(p)
upper_tri_ind <- upper.tri(L)
upper_tri_ind2 <- upper.tri(L, TRUE)
colmn_ind <- col(upper_tri_ind)[upper_tri_ind]
diags <- cbind(1:p, 1:p)
diags2 <- cbind(2:p, 2:p)

b <- MASS::mvrnorm(K * 20, rep(0, p), D)

M <- 3000L
acceptance_sds <- res_sds <- matrix(0.0, M, p)
scale_sds <- rep(0.1, p)
acceptance_L <- matrix(0.0, M, K)
res_L <- matrix(0.0, M, K + p)
scale_L <- rep(0.1, K)
#
current_sds <- sds
current_L <- L

system.time({
    for (m in seq_len(M)) {
        # update the sds
        t_inv_current_L <- backsolve(r = current_L, x = Ip, transpose = TRUE)
        for (i in seq_len(p)) {
            current_sds_i <- current_sds[i]
            scale_sds_i <- scale_sds[i]
            log_mu_i <- log(current_sds_i) - 0.5 * scale_sds_i^2
            proposed_sds_i <- rlnorm(1L, log_mu_i, scale_sds_i)
            pr_sds <- current_sds
            pr_sds[i] <- proposed_sds_i
            numerator_sds_i <- logPC_D_sds(pr_sds, t_inv_current_L,
                                           half_t_mean = 10 * sds) +
                dlnorm(current_sds_i, log_mu_i, scale_sds_i, log = TRUE)
            if (i == 1) {
                denominator_sds_i <- logPC_D_sds(current_sds, t_inv_current_L,
                                                 half_t_mean = 10 * sds) +
                    dlnorm(proposed_sds_i, log_mu_i, scale_sds_i, log = TRUE)
            }
            log_ratio_i <- numerator_sds_i - denominator_sds_i
            if (log_ratio_i > log(runif(1))) {
                current_sds <- pr_sds
                denominator_sds_i <- numerator_sds_i
                acceptance_sds[m, i] <- 1
            }
            res_sds[m, i] <- current_sds[i]
            if (m > 31) {
                scale_sds[i] <- robbins_monro_univ(scale = scale_sds_i,
                                                   acceptance_it = acceptance_sds[m, i],
                                                   it = m - 30)
            }
        }
        # update the off-diagonal elements of the L matrix
        for (i in seq_len(K)) {
            current_L_i <- current_L[upper_tri_ind][i]
            scale_L_i <- scale_L[i]
            proposed_L_i <- runif(1L, min = current_L_i - 0.5 * scale_L_i * sqrt(12),
                                  max = current_L_i + 0.5 * scale_L_i * sqrt(12))
            pr_L <- current_L
            pr_L[upper_tri_ind][i] <- proposed_L_i
            # to ensure that L is the Cholesky factor of the correlation matrix R
            # we need that the diagonal elements are positive, and that the columns
            # have a unit Euclidean length. That is, if x are the elements of one
            # column, then sqrt(sum(x * x)) should be 1. Below we define the diagonal
            # elements as such. But it can happen that the proposed_L_i does not satisfy
            # this contraint, i.e., that sum(x * x) > 1. In this case, we should not
            # accept this L matrix
            ll <- pr_L[seq(1, colmn_ind[i] - 1), colmn_ind[i]]
            ss <- sum(ll * ll)
            if (ss > 1) {
                numerator_L_i <- -1e15
            } else {
                pr_L[colmn_ind[i], colmn_ind[i]] <- sqrt(1 - ss)
                numerator_L_i <- logPC_D_L(pr_L, current_sds)
            }
            if (i == 1) denominator_L_i <- logPC_D_L(current_L, current_sds)
            log_ratio_i <- numerator_L_i - denominator_L_i
            if (is.finite(log_ratio_i) && log_ratio_i > log(runif(1))) {
                current_L <- pr_L
                denominator_L_i <- numerator_L_i
                acceptance_L[m, i] <- 1
            }
            if (m > 31) {
                scale_L[i] <- robbins_monro_univ(scale = scale_L_i,
                                                 acceptance_it = acceptance_L[m, i],
                                                 it = m - 30)
            }
        }
        res_L[m, ] <- current_L[upper_tri_ind2]
    }
})

ar_sds <- colMeans(acceptance_sds[-seq_len(1000L), ])
ar_sds
res_sds <- res_sds[-seq_len(1000L), ]
for (k in seq_len(p)) {
    plot(res_sds[, k], type = "l", ylab = paste("SDS:", k),
         main = ar_sds[k])
}

ar_L <- colMeans(acceptance_L[-seq_len(1000L), ])
ar_L

res_L <- res_L[-seq_len(1000L), ]
for (k in seq_len(K + p)[-1L]) {
    plot(res_L[, k], type = "l", ylab = paste("Element:", k),
         main = ar_L[k])
}

mean_D <- mapply(reconstr_D, L = split(res_L, row(res_L)),
                 sds = split(res_sds, row(res_sds)), SIMPLIFY = FALSE)
mean_D <- Reduce("+", mean_D) / length(mean_D)

# check posterior mean
round(cbind(mean_D[upper.tri(mean_D, TRUE)], D[upper.tri(D, TRUE)]), 3)


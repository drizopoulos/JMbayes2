library("microbenchmark")
source(file.path(getwd(), "Development/MCMC/Surv_Model/prepare_test_fast_log_surv.R"))


log_density_surv2 <- function (bs_gammas, gammas, alphas) {
    lambda_H <- W0_H %*% bs_gammas + W_H %*% gammas
    for (i in seq_along(Wlong_H)) {
        lambda_H <- lambda_H + Wlong_H[[i]] %*% alphas[[i]]
    }
    lambda_h <- matrix(0.0, n, 1)
    if (length(which_event)) {
        lambda_h <- W0_h %*% bs_gammas + W_h %*% gammas
        for (i in seq_along(Wlong_h)) {
            W_h_i <- Wlong_h[[i]]
            lambda_h <- lambda_h + W_h_i %*% alphas[[i]]
        }
    }
    lambda_H2 <- matrix(0.0, nrow(Wlong_H2[[1]]), 1)
    if (length(which_interval)) {
        lambda_H2 <- W0_H2 %*% bs_gammas + W_H2 %*% gammas
        for (i in seq_along(Wlong_H2)) {
            W_H2_i <- Wlong_H2[[i]]
            lambda_H2 <- lambda_H2 + W_H2_i %*% alphas[[i]]
        }
    }
    H <- rowsum(exp(log_Pwk + lambda_H), group = id_H[[1]], reorder = FALSE)
    log_Lik_surv <- numeric(n)
    which_right_event <- c(which_right, which_event)
    if (length(which_right_event)) {
        log_Lik_surv[which_right_event] <- - H[which_right_event]
    }
    if (length(which_event)) {
        log_Lik_surv[which_event] <- log_Lik_surv[which_event] + lambda_h[which_event]
    }
    if (length(which_left)) {
        log_Lik_surv[which_left] <- log1p(- exp(- H[which_left]))
    }
    if (length(which_interval)) {
        H2 <- rowsum(exp(log_Pwk2 + lambda_H2), group = id_H2[[1]], reorder = FALSE)
        log_Lik_surv[which_interval] <- log(exp(- H[which_interval]) -
                                                exp(- (H2[which_interval] + H[which_interval])))
    }
    sum(log_Lik_surv, na.rm = TRUE)
}

################################################################################
################################################################################

test1 <- log_density_surv(bs_gammas, gammas, alphas)
test2 <- log_density_surv2(bs_gammas, gammas, alphas)

all.equal(test1, test2)


microbenchmark(
    old = log_density_surv(bs_gammas, gammas, alphas),
    new = log_density_surv2(bs_gammas, gammas, alphas),
    times = 1000
)









library("microbenchmark")
source(file.path(getwd(), "Development/MCMC/Surv_Model/prepare_test_fast_log_surv.R"))

Wlong_alphas_fun <- function (Wlong, alphas) {
    out <- numeric(nrow(Wlong[[1L]]))
    for (i in seq_along(Wlong)) {
        out <- out + Wlong[[i]] %*% alphas[[i]]
    }
    out
}

group_sum <- function (x, ind) {
    xx <- c(0, cumsum(x)[ind])
    xx[-1L] - xx[-length(xx)]
}

log_density_surv2 <- function (W0H_bs_gammas, WH_gammas, WlongH_alphas,
                               W0h_bs_gammas, Wh_gammas, Wlongh_alphas,
                               W0H2_bs_gammas, WH2_gammas, WlongH2_alphas) {
    lambda_H <- W0H_bs_gammas + WH_gammas + WlongH_alphas
    lambda_h <- matrix(0.0, n, 1)
    if (length(which_event)) {
        lambda_h <- W0h_bs_gammas + Wh_gammas + Wlongh_alphas
    }
    lambda_H2 <- matrix(0.0, nrow(Wlong_H2[[1]]), 1)
    if (length(which_interval)) {
        lambda_H2 <- W0H2_bs_gammas + WH2_gammas + WlongH2_alphas
    }
    H <- group_sum(exp(log_Pwk + lambda_H), indFast_H)
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
        H2 <- group_sum(exp(log_Pwk2 + lambda_H2), indFast_H2)
        log_Lik_surv[which_interval] <- log(exp(- H[which_interval]) -
                                exp(- (H2[which_interval] + H[which_interval])))
    }
    sum(log_Lik_surv, na.rm = TRUE)
}

################################################################################
################################################################################

indFast_H <- id_H[[1]]
indFast_H <- c(indFast_H[-length(indFast_H)] != indFast_H[-1L], TRUE)
indFast_H2 <- id_H[[1]]
indFast_H2 <- c(indFast_H2[-length(indFast_H2)] != indFast_H2[-1L], TRUE)
W0H_bs_gammas <- W0_H %*% bs_gammas
WH_gammas <- W_H %*% gammas
WlongH_alphas <- Wlong_alphas_fun(Wlong_H, alphas)
if (length(which_event)) {
    W0h_bs_gammas <- W0_h %*% bs_gammas
    Wh_gammas <- W_h %*% gammas
    Wlongh_alphas <- Wlong_alphas_fun(Wlong_h, alphas)
}
if (length(which_interval)) {
    W0H2_bs_gammas <- W0_H2 %*% bs_gammas
    WH2_gammas <- W_H2 %*% gammas
    WlongH2_alphas <- Wlong_alphas_fun(Wlong_H2, alphas)
}

test1 <- log_density_surv(bs_gammas, gammas, alphas)
test2 <- log_density_surv2(W0H_bs_gammas, WH_gammas, WlongH_alphas,
                           W0h_bs_gammas, Wh_gammas, Wlongh_alphas,
                           W0H2_bs_gammas, WH2_gammas, WlongH2_alphas)

all.equal(test1, test2)


microbenchmark(
    old = log_density_surv(bs_gammas, gammas, alphas),
    new = {
        W_H %*% gammas
        if (length(which_event)) Wh_gammas <- W_h %*% gammas
        if (length(which_interval)) WH2_gammas <- W_H2 %*% gammas
        log_density_surv2(W0H_bs_gammas, WH_gammas, WlongH_alphas,
                            W0h_bs_gammas, Wh_gammas, Wlongh_alphas,
                            W0H2_bs_gammas, WH2_gammas, WlongH2_alphas)
        },
    times = 2000
)









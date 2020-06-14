source(file.path(getwd(), "Development/MCMC/Surv_Model/examples.R"))
source(file.path(getwd(), "Development/MCMC/Surv_Model/sample_Surv_Funs.R"))

M <- 5000L
res_bs_gammas <- acceptance_bs_gammas <- matrix(0.0, M, length(bs_gammas))
vcov_prop_bs_gammas <- test$vcov_prop$vcov_prop_bs_gammas
scale_bs_gammas <- rep(0.1, length(bs_gammas))
prior_mean_bs_gammas <- test$priors$mean_bs_gammas
prior_Tau_bs_gammas <- test$priors$Tau_bs_gammas
post_A_tau_bs_gammas <- test$priors$A_tau_bs_gammas +
    0.5 * test$priors$rank_Tau_bs_gammas
prior_B_tau_bs_gammas <- test$priors$B_tau_bs_gammas
res_tau_bs_gammas <- numeric(M)
#
any_gammas <- test$model_data$any_gammas
res_gammas <- acceptance_gammas <- matrix(0.0, M, length(gammas))
vcov_prop_gammas <- test$vcov_prop$vcov_prop_gammas
scale_gammas <- rep(0.1, length(gammas))
prior_mean_gammas <- test$priors$mean_gammas
#
res_alphas <- acceptance_alphas <- lapply(alphas,
                                          function (a) matrix(0, M, length(a)))
vcov_prop_alphas <- test$vcov_prop$vcov_prop_alphas
scale_alphas <- lapply(alphas, function (a) a * 0 + 0.1)
prior_mean_alphas <- unlist(test$priors$mean_alphas, use.names = FALSE)
####
tau_bs_gammas <- 2
current_bs_gammas <- jitter(bs_gammas, 80)
current_gammas <- gammas
current_alphas <- alphas

W_h <- scale(W_h, scale = FALSE)
W_H <- scale(W_H, scale = FALSE)
W_H2 <- scale(W_H2, scale = FALSE)

system.time({
    for (m in seq_len(M)) {
        if (m == 1) denominator_surv <- logPC_surv(current_bs_gammas, current_gammas,
                                                   current_alphas, tau_bs_gammas)
        # Update bs_gammas
        for (i in seq_along(current_bs_gammas)) {
            proposed_bs_gammas <- current_bs_gammas
            proposed_bs_gammas[i] <- rnorm(1L, current_bs_gammas[i],
                                           scale_bs_gammas[i])
            numerator_surv <- logPC_surv(proposed_bs_gammas, current_gammas,
                                         current_alphas, tau_bs_gammas)
            log_ratio <- numerator_surv - denominator_surv
            if (is.finite(log_ratio) && min(1, exp(log_ratio)) > runif(1)) {
                current_bs_gammas <- proposed_bs_gammas
                denominator_surv <- numerator_surv
                acceptance_bs_gammas[m, i] <- 1
            }
            if (m > 20) {
                scale_bs_gammas[i] <-
                    robbins_monro_univ(scale = scale_bs_gammas[i],
                                       acceptance_it = acceptance_bs_gammas[m, i],
                                       it = m, target_acceptance = 0.45)
            }
        }
        post_B_tau_bs_gammas <- prior_B_tau_bs_gammas +
            0.5 * c(crossprod(current_bs_gammas, prior_Tau_bs_gammas) %*%
                        current_bs_gammas)
        tau_bs_gammas <- rgamma(1L, post_A_tau_bs_gammas, post_B_tau_bs_gammas)
        # Update gammas
        if (any_gammas) {
            for (i in seq_along(current_gammas)) {
                proposed_gammas <- current_gammas
                proposed_gammas[i] <- rnorm(1L, current_gammas[i],
                                               scale_gammas[i])
                numerator_surv <- logPC_surv(current_bs_gammas, proposed_gammas,
                                             current_alphas, tau_bs_gammas)
                log_ratio <- numerator_surv - denominator_surv
                if (is.finite(log_ratio) && min(1, exp(log_ratio)) > runif(1)) {
                    current_gammas <- proposed_gammas
                    denominator_surv <- numerator_surv
                    acceptance_gammas[m, i] <- 1
                }
                if (m > 20) {
                    scale_gammas[i] <-
                        robbins_monro_univ(scale = scale_gammas[i],
                                           acceptance_it = acceptance_gammas[m, i],
                                           it = m, target_acceptance = 0.45)
                }
            }
        }
        ###
        # updates alphas
        for (i in seq_along(current_alphas)) {
            for (j in seq_along(current_alphas[[i]])) {
                proposed_alphas <- current_alphas
                proposed_alphas[[i]][j] <- rnorm(1L, current_alphas[[i]][j],
                                                 scale_alphas[[i]][j])
                numerator_surv <- logPC_surv(current_bs_gammas, current_gammas,
                                             proposed_alphas, tau_bs_gammas)
                log_ratio <- numerator_surv - denominator_surv
                if (is.finite(log_ratio) && min(1, exp(log_ratio)) > runif(1)) {
                    current_alphas <- proposed_alphas
                    denominator_surv <- numerator_surv
                    acceptance_alphas[[i]][m, j] <- 1
                }
                if (m > 20) {
                    scale_alphas[i] <-
                        robbins_monro_univ(scale = scale_alphas[[i]][j],
                                           acceptance_it = acceptance_alphas[[i]][m, j],
                                           it = m, target_acceptance = 0.45)
                }
                res_alphas[[i]][m, j] <- current_alphas[[i]][j]
            }
        }
        ###
        res_bs_gammas[m, ] <- current_bs_gammas
        res_tau_bs_gammas[m] <- tau_bs_gammas
        res_gammas[m, ] <- current_gammas
        ###
    }
})

###########################

ar_bs_gammas <- colMeans(acceptance_bs_gammas[-seq_len(1000L), ])
ar_bs_gammas

ar_gammas <- colMeans(acceptance_gammas[-seq_len(1000L), , drop = FALSE])
ar_gammas

res_bs_gammas <- res_bs_gammas[-seq_len(1000L), ]
for (k in seq_len(length(current_bs_gammas))) {
    plot(res_bs_gammas[, k], type = "l", ylab = paste0("bs_gammas", k))
}

plot(res_tau_bs_gammas[-seq_len(1000L)], type = "l")

res_gammas <- res_gammas[-seq_len(1000L), , drop = FALSE]
for (k in seq_len(length(current_gammas))) {
    plot(res_gammas[, k], type = "l", ylab = paste0("gammas", k))
}


#out1 <- res_bs_gammas
#out2 <- res_bs_gammas
out3 <- res_bs_gammas

for (k in seq_len(length(current_bs_gammas))) {
    v <- cbind(out1[, k], out2[, k], out3[, k])
    matplot(v, type = "l", ylab = paste0("bs_gammas", k),
            lty = 1)
}





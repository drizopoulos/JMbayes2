source(file.path(getwd(), "Development/MCMC/Surv_Model/examples.R"))
source(file.path(getwd(), "Development/MCMC/Surv_Model/sample_Surv_Funs.R"))

M <- 5000L
res_bs_gammas <- acceptance_bs_gammas <- matrix(0.0, M, length(bs_gammas))
vcov_prop_bs_gammas <- test$vcov_prop$vcov_prop_bs_gammas
scale_bs_gammas <- rep(0.1, length(bs_gammas))
#
prior_mean_bs_gammas <- test$priors$mean_bs_gammas
prior_Tau_bs_gammas <- test$priors$Tau_bs_gammas
post_A_tau_bs_gammas <- test$priors$A_tau_bs_gammas +
    0.5 * test$priors$rank_Tau_bs_gammas
prior_B_tau_bs_gammas <- test$priors$B_tau_bs_gammas
res_tau_bs_gammas <- numeric(M)
####
tau_bs_gammas <- 2
current_bs_gammas <- jitter(bs_gammas, 80)
current_gammas <- gammas
current_alphas <- alphas

system.time({
    for (m in seq_len(M)) {
        if (m == 1) denominator_surv <- logPC_surv(current_bs_gammas, current_gammas,
                                                   current_alphas, tau_bs_gammas)
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
        ###
        res_bs_gammas[m, ] <- current_bs_gammas
        res_tau_bs_gammas[m] <- tau_bs_gammas
        ###
    }
})

###########################

ar_bs_gammas <- mean(acceptance_bs_gammas[-seq_len(1000L)])
ar_bs_gammas

res_bs_gammas <- res_bs_gammas[-seq_len(1000L), ]
for (k in seq_len(length(current_bs_gammas))) {
    plot(res_bs_gammas[, k], type = "l", ylab = paste0("bs_gammas", k))
}

plot(res_tau_bs_gammas[-seq_len(1000L)], type = "l")


#out1 <- res_bs_gammas
#out2 <- res_bs_gammas
out3 <- res_bs_gammas

for (k in seq_len(length(current_bs_gammas))) {
    v <- cbind(out1[, k], out2[, k], out3[, k])
    matplot(v, type = "l", ylab = paste0("bs_gammas", k),
            lty = 1)
}





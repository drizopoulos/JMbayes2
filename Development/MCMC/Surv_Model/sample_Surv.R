source(file.path(getwd(), "Development/MCMC/Surv_Model/examples.R"))
source(file.path(getwd(), "Development/MCMC/Surv_Model/sample_Surv_Funs.R"))

log_density_surv(bs_gammas + 0.1, gammas, alphas)


M <- 5000L
acceptance_bs_gammas <- numeric(M)
res_bs_gammas <- matrix(0.0, M, length(bs_gammas))
vcov_prop_bs_gammas <- test$vcov_prop$vcov_prop_bs_gammas
scale_bs_gammas <- 0.1
#
current_bs_gammas <- bs_gammas
current_gammas <- gammas
current_alphas <- alphas


for (m in seq_len(M)) {
    if (i == 1) denominator_surv <- log_density_surv(current_bs_gammas,
                                                     current_gammas,
                                                     current_alphas)
    proposed_bs_gammas <- rmvnorm(1L, mu = current_bs_gammas,
                                  Sigma = scale_bs_gammas * vcov_prop_bs_gammas)
    numerator_surv <- log_density_surv(proposed_bs_gammas,
                                       current_gammas,
                                       current_alphas)
}

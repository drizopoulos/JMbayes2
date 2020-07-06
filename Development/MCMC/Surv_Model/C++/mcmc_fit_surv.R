#mcmc_fit <- function (model_data, model_info, initial_values, priors, control) {
    # initial values


id_H <- id_H2 <- rep(seq_len(model_data$n), each = control$GK_k)

model_data <- c(model_data, create_Wlong_mats(model_data, model_info,
                                              initial_values, priors,
                                              control),
                list(id_H = id_H, id_H2 = id_H2))

priors$mean_alphas <- unlist(priors$mean_alphas, use.names = FALSE)
priors$Tau_alphas <- .bdiag(priors$Tau_alphas)

system.time(xxx <- mcmc(model_data, model_info, initial_values, priors, control))


#}

for (k in 0:12) {
    if (k == 0)
        plot(xxx$mcmc$tau_bs_gammas, type = "l")
    else
        plot(xxx$mcmc$bs_gammas[, k], type = "l")
}

#mcmc_fit <- function (model_data, model_info, initial_values, priors, control) {
    # initial values


id_H <- id_H2 <- rep(seq_len(model_data$n), each = control$GK_k)

model_data <- c(model_data, create_Wlong_mats(model_data, model_info,
                                              initial_values, priors,
                                              control),
                list(id_H = id_H, id_H2 = id_H2))

priors$mean_alphas <- unlist(priors$mean_alphas, use.names = FALSE)
priors$Tau_alphas <- .bdiag(priors$Tau_alphas)

xxx <- mcmc(model_data, model_info, initial_values, priors, control)

all.equal((xxx$WlongH_alphas), WlongH_alphas, check.attributes = FALSE)

#}

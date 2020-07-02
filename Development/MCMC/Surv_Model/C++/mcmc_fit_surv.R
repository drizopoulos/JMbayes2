#mcmc_fit <- function (model_data, model_info, initial_values, priors, control) {
    # initial values


id_H <- id_H2 <- rep(seq_len(model_data$n), each = control$GK_k)

model_data <- c(model_data, create_Wlong_mats(model_data, model_info,
                                              initial_values, priors,
                                              control),
                list(id_H = id_H, id_H2 = id_H2))



xxx <- mcmc(model_data, model_info, initial_values, priors, control)

all.equal(c(xxx$id_H), which(indFast_H) - 1, check.attributes = FALSE)

#}

#mcmc_fit <- function (model_data, model_info, initial_values, priors, control) {
    # initial values


id_H <- id_H2 <- rep(seq_len(model_data$n), each = control$GK_k)

model_data <- c(model_data, create_Wlong_mats(model_data, model_info,
                                              initial_values, priors,
                                              control),
                list(id_H = id_H, id_H2 = id_H2))


model_data$W_H <- scale(model_data$W_H, scale = FALSE)
model_data$W_h <- scale(model_data$W_h, scale = FALSE)
model_data$W_H2 <- scale(model_data$W_H2, scale = FALSE)
model_data$W_bar <- rbind(attr(model_data$W_h, "scaled:center"))

model_data$Wlong_H <- lapply(model_data$Wlong_H, scale, scale = FALSE)
model_data$Wlong_h <- lapply(model_data$Wlong_h, scale, scale = FALSE)
model_data$Wlong_H2 <- lapply(model_data$Wlong_H2, scale, scale = FALSE)
model_data$Wlong_bar <- lapply(model_data$Wlong_h,
                               function (w) rbind(attr(w, "scaled:center")))


priors$mean_alphas <- unlist(priors$mean_alphas, use.names = FALSE)
priors$Tau_alphas <- .bdiag(priors$Tau_alphas)

system.time(xxx <- mcmc(model_data, model_info, initial_values, priors,
                        control))


#}

for (k in 0:12) {
    if (k == 0)
        plot(xxx$mcmc$tau_bs_gammas, type = "l")
    else
        plot(xxx$mcmc$bs_gammas[, k], type = "l")
}

plot(xxx$mcmc$gammas[, 1], type = "l")
plot(xxx$mcmc$gammas[, 2], type = "l")

plot(xxx$mcmc$alphas[[1]][, 1], type = "l")

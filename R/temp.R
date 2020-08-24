# model_data = Data
# control = con

create_Wlong_mats <- function (model_data, model_info, initial_values, priors,
                               control) {
    betas <- initial_values$betas
    b <- initial_values$b
    gammas <- initial_values$gammas
    bs_gammas <- initial_values$bs_gammas
    alphas <- initial_values$alphas
    # outcome vectors and design matrices
    Time_right <- model_data$Time_right
    Time_left <- model_data$Time_left
    Time_start <- model_data$Time_start
    delta <- model_data$delta
    which_event <- model_data$which_event
    which_right <- model_data$which_right
    which_left <- model_data$which_left
    which_interval <- model_data$which_interval
    W0_H <- model_data$W0_H
    W_H <- model_data$W_H
    X_H <- model_data$X_H
    Z_H <- model_data$Z_H
    U_H <- model_data$U_H
    W0_h <- model_data$W0_h
    W_h <- model_data$W_h
    X_h <- model_data$X_h
    Z_h <- model_data$Z_h
    U_h <- model_data$U_h
    W0_H2 <- model_data$W0_H2
    W_H2 <- model_data$W_H2
    X_H2 <- model_data$X_H2
    Z_H2 <- model_data$Z_H2
    U_H2 <- model_data$U_H2
    # other information
    n <- model_data$n
    idT <- model_data$idT
    log_Pwk <- model_data$log_Pwk
    log_Pwk2 <- model_data$log_Pwk2
    functional_forms_per_outcome <-
        model_info$fun_forms$functional_forms_per_outcome
    # id_H is used to repeat the random effects of each subject GK_k times
    id_H <- lapply(X_H, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
    # this is the linear predictor for the longitudinal outcomes evaluated at the
    # Gauss-Kronrod quadrature points
    eta_H <- linpred_surv(X_H, betas, Z_H, b, id_H)
    # Wlong is the design matrix of all longitudinal outcomes according to the specified
    # functional forms per outcome already multiplied with the interaction terms matrix U
    Wlong_H <- create_Wlong(eta_H, functional_forms_per_outcome, U_H)
    if (length(which_event)) {
        id_h <- lapply(X_h, function (x) seq_len(nrow(x[[1]])))
        eta_h <- linpred_surv(X_h, betas, Z_h, b, id_h)
        Wlong_h <- create_Wlong(eta_h, functional_forms_per_outcome, U_h)
    } else {
        Wlong_h <- rep(list(matrix(0.0, length(Time_right), 1)), length(Wlong_H))
    }
    if (length(which_interval)) {
        id_H2 <- lapply(X_H2, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
        eta_H2 <- linpred_surv(X_H2, betas, Z_H, b, id_H2)
        Wlong_H2 <- create_Wlong(eta_H2, functional_forms_per_outcome, U_H2)
    } else {
        Wlong_H2 <- rep(list(matrix(0.0, length(Time_right), 1)), length(Wlong_H))
    }
    list(Wlong_H = Wlong_H, Wlong_h = Wlong_h, Wlong_H2 = Wlong_H2)
}

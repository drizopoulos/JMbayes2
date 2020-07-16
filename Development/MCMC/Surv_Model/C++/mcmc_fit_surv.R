mcmc_fit <- function (model_data, model_info, initial_values, priors, control) {
    id_H <- id_H2 <- rep(seq_len(model_data$n), each = control$GK_k)
    model_data <- c(model_data, create_Wlong_mats(model_data, model_info,
                                                  initial_values, priors,
                                                  control),
                    list(id_H = id_H, id_H2 = id_H2))

    # center the design matrices for the baseline covariates and
    # the longitudinal process
    model_data$W_H <- scale(model_data$W_H, scale = FALSE)
    model_data$W_h <- scale(model_data$W_h, scale = FALSE)
    model_data$W_H2 <- scale(model_data$W_H2, scale = FALSE)
    model_data$W_bar <- rbind(attr(model_data$W_h, "scaled:center"))

    model_data$Wlong_H <- lapply(model_data$Wlong_H, scale, scale = FALSE)
    model_data$Wlong_h <- lapply(model_data$Wlong_h, scale, scale = FALSE)
    model_data$Wlong_H2 <- lapply(model_data$Wlong_H2, scale, scale = FALSE)
    model_data$Wlong_bar <- lapply(model_data$Wlong_h,
                                   function (w) rbind(attr(w, "scaled:center")))
    # unlist priors and initial values for alphas
    initial_values$alphas <- unlist(initial_values$alphas, use.names = FALSE)
    priors$mean_alphas <- unlist(priors$mean_alphas, use.names = FALSE)
    priors$Tau_alphas <- .bdiag(priors$Tau_alphas)
    # random seed
    if (!exists(".Random.seed", envir = .GlobalEnv))
        runif(1)
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
    n_chains <- control$n_chains
    #n_chains <- 1
    if (n_chains > 1) {
        mcmc_parallel <- function (chain, model_data, model_info, initial_values,
                                   priors, control) {
            seed_ <- control$seed + chain
            set.seed(seed_)
            initial_values[] <- lapply(initial_values, jitter2, factor = 40)
            mcmc_cpp(model_data, model_info, initial_values, priors, control)
        }
        cores <- control$cores
        chains <- split(seq_len(n_chains),
                        rep(seq_len(cores), each = ceiling(n_chains / cores),
                            length.out = n_chains))
        cores <- min(cores, length(chains))
        #cl <- parallel::makeCluster(cores)
        #out <- parallel::parLapply(cl, chains, mcmc_parallel,
        #                           model_data = model_data, model_info = model_info,
        #                           initial_values = initial_values,
        #                           priors = priors, control = control)
        #parallel::stopCluster(cl)
        out <- lapply(chains, mcmc_parallel, model_data = model_data,
                      model_info = model_info, initial_values = initial_values,
                      priors = priors, control = control)
    } else {
        set.seed(control$seed)
        out <- list(mcmc_cpp(model_data, model_info, initial_values, priors,
                             control))
    }
    convert2_mcmclist <- function (name) {
        as.mcmc.list(lapply(out, function (x) as.mcmc(x$mcmc[[name]])))
    }
    list(
        "mcmc" = list(
            "bs_gammas" = convert2_mcmclist("bs_gammas"),
            "tau_bs_gammas" = convert2_mcmclist("tau_bs_gammas"),
            "gammas" = convert2_mcmclist("gammas"),
            "alphas" = convert2_mcmclist("alphas")
        )
    )
}

xxx <- mcmc_fit(model_data, model_info, initial_values, priors, control)

parm <- xxx$mcmc$gammas

gelman.diag(parm)
summary(parm)
effectiveSize(parm)
raftery.diag(parm)
traceplot(parm)
densityplot(parm)
autocorr.plot(parm)

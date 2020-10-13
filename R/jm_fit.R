jm_fit <- function (model_data, model_info, initial_values, priors, control, vcov_prop) {
    # extract family names
    model_info$family_names <- sapply(model_info$families, "[[", "family")
    # extract link names
    model_info$links <- sapply(model_info$families, "[[", "link")
    # set each y element to a matrix
    model_data$y[] <- lapply(model_data$y, as.matrix)
    # for family = binomial and when y has two columns, set the second column
    # to the number of trials instead the number of failures
    binomial_data <- model_info$family_names == "binomial"
    trials_fun <- function (y) {
        if (NCOL(y) == 2) y[, 2] <- y[, 1] + y[, 2]
        y
    }
    model_data$y[binomial_data] <-
        lapply(model_data$y[binomial_data], trials_fun)
    idT <- model_data$idT
    id_H <- rep(paste0(idT, "_", model_data$strata), each = control$GK_k)
    id_H <- match(id_H, unique(id_H))
    id_H_ <- rep(idT, each = control$GK_k)
    id_H_ <- match(id_H_, unique(id_H_))
    id_h <- unclass(idT)
    model_data <-
        c(model_data, create_Wlong_mats(model_data, model_info,
                                        initial_values, priors,
                                        control),
          list(id_H = id_H, id_H_ = id_H_, id_h = id_h))
    # cbind the elements of X_H and Z_H, etc.
    model_data$X_H[] <- lapply(model_data$X_H, docall_cbind)
    model_data$X_h[] <- lapply(model_data$X_h, docall_cbind)
    model_data$X_H2[] <- lapply(model_data$X_H2, docall_cbind)
    model_data$Z_H[] <- lapply(model_data$Z_H, docall_cbind)
    model_data$Z_h[] <- lapply(model_data$Z_h, docall_cbind)
    model_data$Z_H2[] <- lapply(model_data$Z_H2, docall_cbind)
    # center the design matrices for the baseline covariates and
    # the longitudinal process
    model_data$W_bar <- rbind(colMeans(model_data$W_H))
    model_data$W_H <- center_fun(model_data$W_H, model_data$W_bar)
    model_data$W_h <- center_fun(model_data$W_h, model_data$W_bar)
    model_data$W_H2 <- center_fun(model_data$W_H2, model_data$W_bar)

    model_data$Wlong_bar <- lapply(model_data$Wlong_H, colMeans)
    model_data$Wlong_H <- mapply2(center_fun, model_data$Wlong_H,
                                 model_data$Wlong_bar)
    model_data$Wlong_h <- mapply2(center_fun, model_data$Wlong_h,
                                 model_data$Wlong_bar)
    model_data$Wlong_H2 <- mapply2(center_fun, model_data$Wlong_H2,
                                  model_data$Wlong_bar)
    model_data$Wlong_bar <- lapply(model_data$Wlong_bar, rbind)
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
    tik <- proc.time()
    if (n_chains > 1) {
        mcmc_parallel <- function (chain, model_data, model_info, initial_values,
                                   priors, control, vcov_prop) {
            seed_ <- control$seed + chain
            set.seed(seed_)
            not_D <- !names(initial_values) %in% c("betas", "D", "b")
            initial_values[not_D] <- lapply(initial_values[not_D], jitter2)
            mcmc_cpp(model_data, model_info, initial_values, priors, control,
                     vcov_prop)
        }
        cores <- control$cores
        chains <- split(seq_len(n_chains),
                        rep(seq_len(cores), each = ceiling(n_chains / cores),
                            length.out = n_chains))
        cores <- min(cores, length(chains))
        cl <- parallel::makeCluster(cores)
        out <- parallel::parLapply(cl, chains, mcmc_parallel,
                                   model_data = model_data, model_info = model_info,
                                   initial_values = initial_values,
                                   priors = priors, control = control, vcov_prop = vcov_prop)
        parallel::stopCluster(cl)
        #out <- lapply(chains, mcmc_parallel, model_data = model_data,
        #              model_info = model_info, initial_values = initial_values,
        #              priors = priors, control = control, vcov_prop = vcov_prop)
    } else {
        set.seed(control$seed)
        out <- list(mcmc_cpp(model_data, model_info, initial_values, priors,
                             control, vcov_prop))
    }
    tok <- proc.time()
    # create dummy betas
    if (is.null(out[[1]][["mcmc"]][["betas"]])) {
        for (i in seq_along(out)) {
            M <- nrow(out[[i]][["mcmc"]][["bs_gammas"]])
            K <- length(model_data$X)
            for (j in seq_len(K)) {
                k <- ncol(model_data$X[[j]])
                nmn <- paste0("betas", j)
                out[[i]][["mcmc"]][[nmn]] <- matrix(rnorm(M * k), M, k)
            }
        }
    }
    # reconstruct D matrix
    get_D <- function (x) {
        mapply(reconstr_D, split(x$L, row(x$L)), split(x$sds, row(x$sds)),
               SIMPLIFY = FALSE)
    }
    for (i in seq_along(out)) {
        out[[i]][["mcmc"]][["D"]] <-
            do.call("rbind", lapply(get_D(out[[i]][["mcmc"]]), c))
    }
    # Set names
    for (i in seq_along(out)) {
        colnames(out[[i]][["mcmc"]][["bs_gammas"]]) <-
            paste0("bs_gammas_",
                   seq_along(out[[i]][["mcmc"]][["bs_gammas"]][1, ]))
        colnames(out[[i]][["mcmc"]][["tau_bs_gammas"]]) <-
            if ((n_taus <- ncol(out[[i]][["mcmc"]][["tau_bs_gammas"]])) > 1) {
                paste0("tau_bs_gammas_", seq_len(n_taus))
            } else "tau_bs_gammas"
        colnames(out[[i]][["mcmc"]][["gammas"]]) <- colnames(model_data$W_H)
        colnames(out[[i]][["mcmc"]][["alphas"]]) <-
            unlist(lapply(model_data$U_H, colnames), use.names = FALSE)
        ind <- lower.tri(initial_values$D, TRUE)
        colnames(out[[i]][["mcmc"]][["D"]]) <-
            paste0("D[", row(initial_values$D)[ind], ", ",
                   col(initial_values$D)[ind], "]")
        for (j in seq_along(model_data$X)) {
            colnames(out[[i]][["mcmc"]][[paste0("betas", j)]]) <-
                colnames(model_data$X[[j]])
        }
        colnames(out[[i]][["mcmc"]][["sigmas"]]) <-
            paste0("sigmas_",
                   seq_along(out[[i]][["mcmc"]][["sigmas"]][1, ]))
    }
    # drop sigmas that are not needed
    has_sigmas <- initial_values$log_sigmas > -20.0
    for (i in seq_along(out)) {
        out[[i]][["mcmc"]][["sigmas"]] <-
            out[[i]][["mcmc"]][["sigmas"]][, has_sigmas, drop = FALSE]
    }
    convert2_mcmclist <- function (name) {
        as.mcmc.list(lapply(out, function (x) {
            kk <- x$mcmc[[name]]
            if (length(d <- dim(kk)) > 2) {
                m <- matrix(0.0, d[3L], d[1L] * d[2L])
                for (j in seq_len(d[3L])) m[j, ] <- c(kk[, , j])
                as.mcmc(m)
            } else as.mcmc(kk)
        }))
    }
    get_acc_rates <- function (name_parm) {
        do.call("rbind", lapply(out, function (x) x[["acc_rate"]][[name_parm]]))
    }
    if (control$save_random_effects) {
        parms <- c("bs_gammas", "tau_bs_gammas", "gammas", "alphas", "W_bar_gammas",
                   "Wlong_bar_alphas", "D", paste0("betas", seq_along(model_data$X)),
                   "sigmas", "b")
    } else {
        parms <- c("bs_gammas", "tau_bs_gammas", "gammas", "alphas", "W_bar_gammas",
                   "Wlong_bar_alphas", "D", paste0("betas", seq_along(model_data$X)),
                   "sigmas")
    }
    if (!length(attr(model_info$terms$terms_Surv_noResp, "term.labels")))
        parms <- parms[parms != "gammas"]
    mcmc_out <- lapply_nams(parms, convert2_mcmclist)
    if (control$save_random_effects) {
        list(
            "mcmc" = mcmc_out,
            "acc_rates" = lapply_nams(parms, get_acc_rates),
            "running_time" = tok - tik
        )
    } else {
        postmeans_b = Reduce('+', lapply(out, function(x) x$mcmc$b[, , 1])) / n_chains
        postvars_b = Reduce('+', lapply(out, function(x) x$mcmc$var_b)) / n_chains
        list(
            "postmeans_b" = postmeans_b,
            'postvars_b' = postvars_b,
            "mcmc" = mcmc_out,
            "acc_rates" = lapply_nams(parms, get_acc_rates),
            "running_time" = tok - tik
        )
    }
}

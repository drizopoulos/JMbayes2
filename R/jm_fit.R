jm_fit <- function (model_data, model_info, initial_values, priors, control) {
    # extract family names
    model_info$family_names <- sapply(model_info$families, "[[", "family")
    # extract link names
    model_info$links <- sapply(model_info$families, "[[", "link")
    # set each y element to a matrix
    model_data$y[] <- lapply(model_data$y, as.matrix)
    # for family = binomial and when y has two columns, set the second column
    # to the number of trials instead the number of failures
    binomial_data <- model_info$family_names %in% c("binomial", "beta binomial")
    trials_fun <- function (y) {
        if (NCOL(y) == 2L) y[, 2L] <- y[, 1L] + y[, 2L]
        y
    }
    model_data$y[binomial_data] <-
        lapply(model_data$y[binomial_data], trials_fun)
    idT <- model_data$idT
    # When we have competing risks we have the same idT for the different strata
    # id_H below distinguishes between the different risks per subject. Then for
    # each unique idT & strata/risk combination we have the quadrature points.
    ni_event <- tapply(idT, idT, length)
    model_data$ni_event <- cbind(c(0, head(cumsum(ni_event), -1)),
                                 cumsum(ni_event))
    id_H <- rep(paste0(idT, "_", unlist(tapply(idT, idT, seq_along))),
                each = control$GK_k)
    id_H <- match(id_H, unique(id_H))
    # id_H_ repeats each unique idT the number of quadrature points
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
    # center design matrix fixed effects
    model_data$X[] <- mapply2(center_fun, model_data$X, model_data$Xbar)
    model_data$Xbar <- lapply(model_data$Xbar, rbind)
    # center the design matrices for the baseline covariates and
    # the longitudinal process
    colSds <- function (m) apply(m, 2L, sd)
    model_data$W_bar <- rbind(colMeans(model_data$W_H))
    model_data$W_sds <- rbind(colSds(model_data$W_H))
    model_data$W_std <- model_data$W_bar / model_data$W_sds
    model_data$W_H <- center_fun(model_data$W_H, model_data$W_bar,
                                 model_data$W_sds)
    model_data$W_h <- center_fun(model_data$W_h, model_data$W_bar,
                                 model_data$W_sds)
    model_data$W_H2 <- center_fun(model_data$W_H2, model_data$W_bar,
                                  model_data$W_sds)

    model_data$Wlong_bar <- lapply(model_data$Wlong_H, colMeans)
    model_data$Wlong_sds <- lapply(model_data$Wlong_H, colSds)
    model_data$Wlong_H <- mapply2(center_fun, model_data$Wlong_H,
                                 model_data$Wlong_bar, model_data$Wlong_sds)
    model_data$Wlong_h <- mapply2(center_fun, model_data$Wlong_h,
                                 model_data$Wlong_bar, model_data$Wlong_sds)
    model_data$Wlong_H2 <- mapply2(center_fun, model_data$Wlong_H2,
                                  model_data$Wlong_bar, model_data$Wlong_sds)
    model_data$Wlong_std <- mapply2("/", model_data$Wlong_bar, model_data$Wlong_sds)
    model_data$Wlong_bar <- lapply(model_data$Wlong_bar, rbind)
    model_data$Wlong_sds <- lapply(model_data$Wlong_sds, rbind)
    model_data$Wlong_std <- lapply(model_data$Wlong_std, rbind)
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
    cores <- control$cores
    chains <- seq_len(n_chains)
    mcmc_parallel <- function (chain, model_data, model_info, initial_values,
                               priors, control) {
      not_D <- !names(initial_values) %in% c("D")
      initial_values[not_D] <- lapply(initial_values[not_D], jitter2)
      mcmc_cpp(model_data, model_info, initial_values, priors, control)
    }
    if (cores > 1L) {
        cores <- min(cores, length(chains))
        cl <- parallel::makeCluster(cores)
        parallel::clusterSetRNGStream(cl = cl, iseed = control$seed)
        out <- parallel::parLapply(cl, chains, mcmc_parallel,
                                   model_data = model_data, model_info = model_info,
                                   initial_values = initial_values,
                                   priors = priors, control = control)
        parallel::stopCluster(cl)
    } else {
        set.seed(control$seed)
        out <- lapply(chains, mcmc_parallel, model_data = model_data,
                    model_info = model_info, initial_values = initial_values,
                    priors = priors, control = control)
    }
    tok <- proc.time()
    # split betas per outcome
    ind_FE <- model_data$ind_FE
    for (i in seq_along(out)) {
        betas <- out[[i]][["mcmc"]][["betas"]]
        for (j in seq_along(ind_FE)) {
            nam_j <- paste0("betas", j)
            out[[i]][["mcmc"]][[nam_j]] <- betas[, ind_FE[[j]], drop = FALSE]
        }
        out[[i]][["mcmc"]][["betas"]] <- NULL
    }
    # reconstruct D matrix
    get_D <- function (x) {
        mapply2(reconstr_D, split(x$L, row(x$L)), split(x$sds, row(x$sds)))
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
    keep_its <- seq(1L, control$n_iter - control$n_burnin, by = control$n_thin)
    convert2_mcmclist <- function (name) {
        as.mcmc.list(lapply(out, function (x) {
            kk <- x$mcmc[[name]]
            if (length(d <- dim(kk)) > 2) {
                m <- matrix(0.0, d[3L], d[1L] * d[2L])
                for (j in seq_len(d[3L])) m[j, ] <- c(kk[, , j])
                as.mcmc(m[keep_its, , drop = FALSE])
            } else as.mcmc(kk[keep_its, , drop = FALSE])
        }))
    }
    get_acc_rates <- function (name_parm) {
        do.call("rbind", lapply(out, function (x) x[["acc_rate"]][[name_parm]]))
    }
    parms <- c("bs_gammas", "tau_bs_gammas", "gammas", "alphas", "W_std_gammas",
               "Wlong_std_alphas", "D", paste0("betas", seq_along(model_data$X)),
               "sigmas",
               "alphaF", "frailty") #!! new
    parms2 <- c("bs_gammas", "gammas", "alphas", "L", "sds", "betas", "b",
                "sigmas",
                "alphaF", "frailty") #!! new
    if (control$save_random_effects) parms <- c(parms, "b")
    if (!length(attr(model_info$terms$terms_Surv_noResp, "term.labels")))
        parms <- parms[parms != "gammas"]
    if (all(!has_sigmas))
        parms <- parms[parms != "sigmas"]
    mcmc_out <- lapply_nams(parms, convert2_mcmclist)
    mcmc_out <- list(
        "mcmc" = mcmc_out,
        "acc_rates" = lapply_nams(parms2, get_acc_rates),
        "logLik" = do.call("rbind", lapply(out, "[[", "logLik")),
        "running_time" = tok - tik
    )
    if (!control$save_random_effects) {
        postmeans_b <- Reduce('+', lapply(out, function(x) x$mcmc$b[, , 1])) / n_chains
        cumsum_b <- Reduce('+', lapply(out, function(x) x$mcmc$cumsum_b))
        outprod_b <- Reduce('+', lapply(out, function(x) x$mcmc$outprod_b))
        K <- (control$n_iter - control$n_burnin) * control$n_chains
        means_b <- cumsum_b / K
        outprod_means_b_cube <- array(0.0, dim = dim(outprod_b))
        for (i in 1:nrow(means_b)) {
            outprod_means_b_cube[, , i] <- means_b[i, ] %o% means_b[i, ]
        }
        postvars_b <- (outprod_b / K - outprod_means_b_cube) * K / (K - 1)
        mcmc_out <- c(mcmc_out, list("postmeans_b" = postmeans_b,
                                     'postvars_b' = postvars_b))
    }
    ########################
    # Calculate Statistics
    S <- lapply(mcmc_out$mcmc, summary)
    statistics <- list(
        Mean = lapply(S, get_statistic, "Mean"),
        Median = lapply(S, get_statistic, "Median"),
        SD = lapply(S, get_statistic, "SD"),
        SE = lapply(S, get_statistic, "Time-series SE"),
        CI_low = lapply(S, get_statistic, "2.5CI"),
        CI_upp = lapply(S, get_statistic, "97.5CI"),
        P = lapply(mcmc_out$mcmc, function (x) apply(do.call("rbind", x), 2L, Ptail)),
        Effective_Size = lapply(mcmc_out$mcmc, function (x)
            apply(do.call("rbind", x), 2L, effective_size))
    )
    if (!is.null(mcmc_out$mcmc[["b"]])) {
        znams <- unlist(lapply(model_data$Z, colnames), use.names = FALSE)
        l <- sapply(model_data$unq_idL, length)
        dnames_b <- list(unlist(model_data$unq_idL[which.max(l)]), znams)
        fix_b <- function (stats) {
            x <- stats$b
            dim(x) <- sapply(dnames_b, length)
            dimnames(x) <- dnames_b
            stats$b <- x
            stats
        }
        statistics[] <- lapply(statistics, fix_b)
        nRE <- ncol(statistics$Mean$b)
        b <- do.call("rbind", mcmc_out$mcmc[["b"]])
        jstart <- 1
        jstop <- length(keep_its)
        for (i in seq_len(control$n_chains)) {
          mcmc_out$mcmc[["b"]][[i]] <- array(0.0, dim = c(length(dnames_b[[1]]), nRE, length(keep_its)))
          for (j in jstart:jstop) {
            mcmc_out$mcmc[["b"]][[i]][, , j - ((i - 1) * length(keep_its))] <- b[j, ]
          }
          jstart <- jstart + length(keep_its)
          jstop <- jstop + length(keep_its)
        }
        nn <- nrow(statistics$Mean[["b"]])
        post_vars <- array(0.0, c(nRE, nRE, nn))
        for (i in seq_len(nn)) {
            post_vars[, , i] <- var(b[, seq(0, nRE - 1) * nn + i, drop = FALSE])
        }
        statistics <- c(statistics, post_vars = list(post_vars))
    }
    if (control$n_chains > 1) {
        no_b <- !names(mcmc_out$mcmc) %in% "b"
        calculate_Rhat <- function (theta) {
            ntheta <- ncol(theta[[1]])
            tt <- try(coda::gelman.diag(theta)$psrf, silent = TRUE)
            if (inherits(tt, "try-error"))
                matrix(as.numeric(NA), ntheta, 2,
                       dimnames = list(NULL, c("Point est.", "Upper C.I.")))
            else tt
        }
        statistics <- c(statistics,
                        Rhat = list(lapply(mcmc_out$mcmc[no_b], calculate_Rhat)))
    }
    if (!control$save_random_effects) {
        statistics$Mean$b <- mcmc_out$postmeans_b
        statistics$post_vars <- mcmc_out$postvars_b
        mcmc_out <- mcmc_out[!names(mcmc_out) %in% c('postmeans_b', 'postvars_b')]
        mcmc_out[['mcmc']][['b']] <- lapply(out, function (x) x$mcmc$b_last)
    }
    # Fit statistics
    thetas <- statistics$Mean
    thetas[["betas"]] <- thetas[grep("^betas", names(thetas))]
    thetas[["D"]] <- nearPD(lowertri2mat(thetas[["D"]]))
    if (is.null(thetas[["gammas"]])) thetas[["gammas"]] <- 0.0
    if (is.null(thetas[["sigmas"]])) thetas[["sigmas"]] <- 0.0
    clogLik_mean_parms <- logLik_jm(thetas, model_data, model_info, control)
    conditional_fit_stats <- fit_stats(mcmc_out$logLik, clogLik_mean_parms)
    #
    res_thetas <- thetas
    res_thetas$bs_gammas <- do.call("rbind", mcmc_out$mcmc$bs_gammas)
    res_thetas$gammas <- if (!is.null(mcmc_out$mcmc$gammas)) {
        do.call("rbind", mcmc_out$mcmc$gammas)
    } else matrix(0.0, nrow(res_thetas$bs_gammas), 1)
    res_thetas$alphas <- do.call("rbind", mcmc_out$mcmc$alphas)
    res_thetas$tau_bs_gammas <- do.call("rbind", mcmc_out$mcmc$tau_bs_gammas)
    D <- do.call("rbind", mcmc_out$mcmc$D)
    res_thetas$D <- array(0.0, c(dim(thetas$D), nrow(res_thetas$bs_gammas)))
    for (i in seq_len(nrow(res_thetas$bs_gammas))) {
        res_thetas$D[, , i] <- lowertri2mat(D[i, ])
    }
    res_thetas$sigmas <- if (!is.null(mcmc_out$mcmc$sigmas)) {
        do.call("rbind", mcmc_out$mcmc$sigmas)
    } else matrix(0.0, nrow(res_thetas$bs_gammas), 1)
    res_thetas[["betas"]] <-
        lapply(mcmc_out$mcmc[grep("^betas", names(mcmc_out$mcmc))], function (m)
            do.call("rbind", m))
    res_thetas$frailty <- do.call("rbind", mcmc_out$mcmc$frailty) #!! new
    res_thetas$alphaF  <- do.call("rbind", mcmc_out$mcmc$alphaF) #!! new
    mcmc_out$mlogLik <- mlogLik_jm(res_thetas, statistics$Mean[["b"]],
                          statistics$post_vars, model_data, model_info, control)
    ind <- names(thetas) %in% c("sigmas", "bs_gammas", "gammas", "alphas",
                                "tau_bs_gammas", 
                                "alphaF", "frailty") #!! new
    thetas[ind] <- lapply(thetas[ind], rbind)
    thetas$betas <- lapply(thetas$betas, rbind)
    dim(thetas[["D"]]) <- c(dim(thetas[["D"]]), 1L)
    mlogLik_mean_parms <-
        c(mlogLik_jm(thetas, statistics$Mean[["b"]], statistics$post_vars,
                     model_data, model_info, control))
    marginal_fit_stats <- fit_stats(mcmc_out$mlogLik, mlogLik_mean_parms)
    mcmc_out$logLik <- mcmc_out$mlogLik <- NULL
    c(mcmc_out, list(statistics = statistics,
                     fit_stats = list(conditional = conditional_fit_stats,
                                      marginal = marginal_fit_stats)),
      list(Wlong_bar = model_data$Wlong_bar,
           Wlong_sds = model_data$Wlong_sds, Wlong_std = model_data$Wlong_std,
           W_bar = model_data$W_bar, W_sds = model_data$W_sds,
           W_std = model_data$W_std))
}

opt_model <- function (models, newdata, t0, parallel = "snow", cores = 1L) {
    if (!all(sapply(models, function (obj) inherits(obj, "jm")))) {
        stop("all objects in 'models' must inherit from class 'jm'.\n")
    }
    object <- models[[1L]]
    id_var <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    Time_var <- object$model_info$var_names$Time_var
    TTime_var <- if (length(Time_var) > 1) Time_var[2L] else Time_var
    event_var <- object$model_info$var_names$event_var
    if (!is.data.frame(newdata)) {
        stop("'newdata' must be a data.frame.\n")
    }
    if (!all((c(id_var, time_var) %in% names(newdata)))) {
        stop("'newdata' must contain the variables '", id_var, "' and '",
             time_var, "'.\n")
    }
    ND <- newdata
    ND[[id_var]] <- match(ND[[id_var]], unique(ND[[id_var]]))
    ND_before <- ND[ND[[time_var]] <= t0, ]
    ND_before[[TTime_var]] <- t0
    ND_before[[event_var]] <- 0
    ND_after <- ND[ND[[time_var]] > t0, ]
    data_after_t0 <- nrow(ND_after)
    ###
    MISE_perid <- function (obs, reps, id, times) {
        trapezoid_rule <- function (f, x) {
            sum(0.5 * diff(x) * (f[-length(x)] + f[-1L]))
        }
        loess.smooth2 <- function (x, y) {
            loess.smooth(x, y, span = 0.75, degree = 1,
                         family = "gaussian", evaluation = 200)
        }
        smooth <- function (x, y) {
            n <- length(x)
            if (n >= 4) {
                loess.smooth2(x, y)
            } else if (n > 1 && n < 4) {
                spline(x, y)
            } else {
                list(x = NA_real_, y = NA_real_)
            }
        }
        gof_fun <- function (y, times, id, type) {
            if (type == "variogram") {
                ls <- smooth(times, y)
                ind <- findInterval(times, ls$x)
                rr <- y - ls$y[ind]
                variogram(rr, times, id)[[1L]]
            } else {
                cbind(times, y)
            }
        }
        Obs_ave <- gof_fun(obs, times, id, "average")
        Obs_vario <- gof_fun(obs, times, id, "variogram")
        G_obs_ave <- smooth(Obs_ave[, 1L], Obs_ave[, 2L])
        F_obs_ave <- mapply(smooth, y = split(Obs_ave[, 2L], id),
                            x = split(Obs_ave[, 1L], id), SIMPLIFY = FALSE)
        G_obs_vario <- smooth(Obs_vario[, 1L], Obs_vario[, 2L])
        ni <- tapply(id, id, length)
        id_long <- rep(unique(id), sapply(ni, function (n) if (n < 2) 1 else
            ncol(combn(n, 2))))
        F_obs_vario <- mapply(smooth, y = split(Obs_vario[, 2L], id_long),
                              x = split(Obs_vario[, 1L], id_long),
                              SIMPLIFY = FALSE)
        mise <- function (obs, rep) {
            trapezoid_rule((obs$y - rep$y)^2, obs$x)
        }
        n <- length(unique(id))
        M <- ncol(reps)
        MISE_mod_ave <- MISE_mod_vario <- numeric(M)
        MISE_ave <- MISE_vario <- matrix(0.0, n, M)
        for (m in seq_len(M)) {
            reps_ave <- gof_fun(reps[, m], times, id, "average")
            reps_vario <- gof_fun(reps[, m], times, id, "variogram")
            G_reps_ave <- smooth(reps_ave[, 1L], reps_ave[, 2L])
            F_reps_ave <-
                mapply(smooth, y = split(reps_ave[, 2L], id),
                       x = split(reps_ave[, 1L], id), SIMPLIFY = FALSE)
            G_reps_vario <- smooth(reps_vario[, 1L], reps_vario[, 2L])
            F_reps_vario <-
                mapply(smooth, y = split(reps_vario[, 2L], id_long),
                       x = split(reps_vario[, 1L], id_long), SIMPLIFY = FALSE)
            MISE_mod_ave[m] <- mise(G_obs_ave, G_reps_ave)
            MISE_ave[, m] <- mapply(mise, obs = F_obs_ave, rep = F_reps_ave)
            MISE_mod_vario[m] <- mise(G_obs_vario, G_reps_vario)
            MISE_vario[, m] <- mapply(mise, obs = F_obs_vario, rep = F_reps_vario)
        }
        list(MISE_ave = MISE_ave, MISE_vario = MISE_vario,
             MISE_mod_ave = MISE_mod_ave, MISE_mod_vario = MISE_mod_vario)
    }
    MISE_model <- function (object) {
        prs <- predict(object, newdata = ND_before, newdata2 = ND_after,
                       return_params_mcmc = TRUE)
        sims <- simulate(object, nsim = 200L, newdata = ND_before,
                          include_outcome = TRUE, random_effects = "mcmc",
                          params_mcmc = prs$newdata$params_mcmc)
        n_outcomes <- length(sims[["outcome"]])
        id <- attr(sims$outcome[[1]], "id")
        MISE_mod_ave <- MISE_mod_vario <- matrix(0.0, 200L, n_outcomes)
        MISEs_ave <- MISEs_vario <- array(0, c(length(unique(id)), 200L, n_outcomes))
        dimnames(MISEs_ave) <- dimnames(MISEs_vario) <-
            list(unique(newdata[[id_var]]), colnames(sims[[1]]),
                 names(sims)[seq_len(n_outcomes)])
        for (m in seq_len(n_outcomes)) {
            outcome <- sims$outcome[[m]]
            mises <- MISE_perid(outcome, sims[[m]], attr(outcome, "id"),
                                attr(outcome, "times"))
            MISE_mod_ave[, m] <- mises$MISE_mod_ave
            MISEs_ave[, , m] <- mises$MISE_ave
            MISE_mod_vario[, m] <- mises$MISE_mod_vario
            MISEs_vario[, , m] <- mises$MISE_vario
        }
        list(MISEs_ave = MISEs_ave, MISEs_vario = MISEs_vario, Preds = prs,
             MISE_mod_ave = MISE_mod_ave, MISE_mod_vario = MISE_mod_vario)
    }
    ###
    if (cores > 1L && length(models) > 1L) {
        have_mc <- have_snow <- FALSE
        if (parallel == "multicore") {
            have_mc <- .Platform$OS.type != "windows"
        } else if (parallel == "snow") {
            have_snow <- TRUE
        }
        if (!have_mc && !have_snow) cores <- 1L
        loadNamespace("parallel")
        if (have_mc) {
            out <- parallel::mclapply(models, MISE_model, mc.cores = cores)
        } else {
            cl <- parallel::makePSOCKcluster(rep("localhost", cores))
            invisible(parallel::clusterEvalQ(cl, library("JMbayes2")))
            out <- parallel::parLapply(cl, models, MISE_model)
            parallel::stopCluster(cl)
        }
    } else {
        out <- lapply(models, MISE_model)
    }
    out <- lapply(out, function (x) {
        x$std_MISEs_ave <- x$MISEs_ave
        x$std_MISEs_vario <- x$MISEs_vario
        x
    })
    n_outcomes <- dim(out[[1]][["MISEs_ave"]])[3L]
    for (m in seq_len(n_outcomes)) {
        vals <- do.call("c", lapply(out, function (x) x$MISEs_ave[, , m]))
        out[] <- lapply(out, function (x) {
            x$std_MISEs_ave[, , m] <- (x$MISEs_ave - mean(vals)) / sd(vals)
            x
        })
        vals <- do.call("c", lapply(out, function (x) x$MISEs_vario[, , m]))
        out[] <- lapply(out, function (x) {
            x$std_MISEs_vario[, , m] <- (x$MISEs_vario - mean(vals)) / sd(vals)
            x
        })
    }
    out[] <- lapply(out, function (x) {
        x$MISE_mod_ave <- colMeans(x$MISE_mod_ave, na.rm = TRUE)
        x$MISEs_ave <- apply(x$MISEs_ave, c(1, 3), mean, na.rm = TRUE)
        x$std_MISEs_ave <- apply(x$std_MISEs_ave, c(1, 3), mean, na.rm = TRUE)
        x$MISE_mod_vario <- colMeans(x$MISE_mod_vario, na.rm = TRUE)
        x$MISEs_vario <- apply(x$MISEs_vario, c(1, 3), mean, na.rm = TRUE)
        x$std_MISEs_vario <- apply(x$std_MISEs_vario, c(1, 3), mean, na.rm = TRUE)
        x
    })
    out
}

#object = fm4
#newdata = testing
#newdata2 = testing2
IndvPred_lme <- function (object, newdata, newdata2) {
    if (!inherits(object, "lme") && !inherits(object, "lmeComponents"))
        stop("Use only with 'lme' or 'lmeComponents' objects.\n")
    if (inherits(object, "lme")) {
        TermsX <- object$terms
        formYz <- formula(object$modelStruct$reStruct[[1]])
        mfZ <- model.frame(terms(formYz), data = object$data)
        TermsZ <- attr(mfZ, "terms")
        id_var <- names(object$modelStruct$reStruct)
        betas <- fixef(object)
        sigma <- object$sigma
        D <- lapply(pdMatrix(object$modelStruct$reStruct), "*", sigma^2)[[1]]
        V <- vcov(object)
    } else {
        TermsX <- object$TermsX
        TermsZ <- object$TermsZ
        id_var <- object$idVar
        betas <- object$betas
        sigma <- object$sigma
        D <- object$D
        V <- object$V
    }
    all_vars <- unique(c(all.vars(TermsX), all.vars(TermsZ)))
    newdata_nomiss <- newdata[complete.cases(newdata[all_vars]), ]
    mfX_new <- model.frame(TermsX, data = newdata_nomiss)
    X_new <- model.matrix(TermsX, mfX_new)
    mfZ_new <- model.frame(TermsZ, data = newdata_nomiss)
    Z_new <- model.matrix(TermsZ, mfZ_new)
    y_new <- model.response(mfX_new, "numeric")
    if (length(id_var) > 1)
        stop("the current version of the function only works with a single grouping variable.\n")
    if (is.null(newdata[[id_var]]))
        stop("subject id variable not in newdata.")
    id_nomiss <- match(newdata_nomiss[[id_var]], unique(newdata_nomiss[[id_var]]))
    n <- length(unique(id_nomiss))
    modes <- matrix(0, n, ncol(Z_new))
    post_vars <- DZtVinv <- vector("list", n)
    for (i in seq_len(n)) {
        id_i <- id_nomiss == i
        X_new_id <- X_new[id_i, , drop = FALSE]
        Z_new_id <- Z_new[id_i, , drop = FALSE]
        Vi_inv <- solve(Z_new_id %*% tcrossprod(D, Z_new_id) +
                            sigma^2 * diag(sum(id_i)))
        DZtVinv[[i]] <- tcrossprod(D, Z_new_id) %*% Vi_inv
        modes[i, ] <- c(DZtVinv[[i]] %*% (y_new[id_i] - X_new_id %*%
                                              betas))
        t1 <- DZtVinv[[i]] %*% Z_new_id %*% D
        t2 <- DZtVinv[[i]] %*% X_new_id %*% V %*%
            crossprod(X_new_id, Vi_inv) %*% Z_new_id %*% D
        post_vars[[i]] <- D - t1 + t2
    }
    fitted_y <- c(X_new %*% betas) +
        rowSums(Z_new * modes[id_nomiss, , drop = FALSE])
    #####
    newdata2_nomiss <- newdata2[complete.cases(newdata2[all_vars]), ]
    mfX_new2 <- model.frame(TermsX, data = newdata2_nomiss)
    X_new2 <- model.matrix(TermsX, mfX_new2)
    mfZ_new2 <- model.frame(TermsZ, data = newdata2_nomiss)
    Z_new2 <- model.matrix(TermsZ, mfZ_new2)
    id_nomiss <- match(newdata2_nomiss[[id_var]], unique(newdata2_nomiss[[id_var]]))
    predicted_y <- c(X_new2 %*% betas) +
        rowSums(Z_new2 * modes[id_nomiss, , drop = FALSE])
    list(fitted_y = fitted_y, predicted_y = predicted_y)
}



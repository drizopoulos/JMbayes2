library("JMbayes2")
create_folds <- function (data, V = 5, id_var = "id", seed = 123L) {
    if (!exists(".Random.seed", envir = .GlobalEnv))
        runif(1L)
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
    set.seed(seed)
    data <- as.data.frame(data)
    ids <- data[[id_var]]
    unq_ids <- unique(ids)
    n <- length(unq_ids)
    splits <- split(seq_len(n), sample(rep(seq_len(V), length.out = n)))
    training <- testing <- vector("list", V)
    for (i in seq_along(training)) {
        ind <- ids %in% unq_ids[splits[[i]]]
        training[[i]] <- data[!ind, ]
        testing[[i]] <- data[ind, ]
    }
    list("training" = training, "testing" = testing)
}


if (FALSE) {
    object = Models[[1]][[4]]
    Tstart = 6
    Thoriz = 8
    eps = 0.001
    cores = max(parallel::detectCores() - 1, 1)
    integrated = FALSE
    type_weights = "model-based"
    newdata = aids

    object = Models
    Tstart = tstr
    Thoriz = thor
    eps = 0.001
    cores = max(parallel::detectCores() - 1, 1)
    integrated = TRUE
    model_weights = Brier_weights
    type_weights = "model-based"
    parallel = "snow"
    newdata = test_data
}

tvEPCE <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
                    eps = 0.001, cores = max(parallel::detectCores() - 1, 1),
                    ...) {
    is_jm <- function (object) inherits(object, "jm")
    if (!is_jm(object)) {
        if (!all(sapply(unlist(object, recursive = FALSE), is_jm)))
            stop("Use only with 'jm' objects.\n")
    }
    if (is.null(Thoriz) && is.null(Dt)) {
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    }
    if (!is.null(Thoriz) && Thoriz <= Tstart) {
        stop("'Thoriz' must be larger than 'Tstart'.")
    }
    if (is.null(Thoriz)) {
        Thoriz <- Tstart + Dt
    }
    Tstart <- Tstart + 1e-06
    epce_fun <- function (qi_u_t, qi_u_t2, tilde_event, eps) {
        - mean(tilde_event * (log(1 - qi_u_t2) - log(eps)) + log(qi_u_t))
    }
    if (!is.data.frame(newdata) && !is.list(newdata) &&
        !all(sapply(newdata, is.data.frame))) {
        stop("'newdata' must be a data.frame or a list of data.frames.\n")
    }
    # if newdata is a list and not a data.frame,
    # Super Learning will be used
    if (!is.data.frame(newdata) && is.list(newdata)) {
        folds <- rep(seq_along(newdata), sapply(newdata, nrow))
        newdata <- do.call("rbind", newdata)
        newdata[["fold_"]] <- folds
    }
    # if Super Learning, object needs to be a list with length the
    # number of folds. In each element of the list, we have a list of fitted
    # models
    obj <- if (is_jm(object)) object else object[[1L]][[1L]]
    id_var <- obj$model_info$var_names$idVar
    time_var <- obj$model_info$var_names$time_var
    Time_var <- obj$model_info$var_names$Time_var
    event_var <- obj$model_info$var_names$event_var
    type_censoring <- object$model_info$type_censoring
    if (obj$model_info$CR_MS) {
        stop("'tvEPCE()' currently only works for right censored data.")
    }
    if (is.null(newdata[[id_var]])) {
        stop("cannot find the '", id_var, "' variable in newdata.", sep = "")
    }
    if (is.null(newdata[[time_var]])) {
        stop("cannot find the '", time_var, "' variable in newdata.", sep = "")
    }
    if (any(sapply(Time_var, function (nmn) is.null(newdata[[nmn]])))) {
        stop("cannot find the '", paste(Time_var, collapse = ", "),
             "' variable(s) in newdata.", sep = "")
    }
    if (is.null(newdata[[event_var]])) {
        stop("cannot find the '", event_var, "' variable in newdata.", sep = "")
    }
    newdata <- newdata[newdata[[Time_var]] > Tstart, ]
    newdata <- newdata[newdata[[time_var]] <= Tstart, ]
    if (!nrow(newdata)) {
        stop("there are no data on subjects who had an observed event time after Tstart ",
             "and longitudinal measurements before Tstart.")
    }
    newdata[[id_var]] <- newdata[[id_var]][, drop = TRUE]
    test <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
    if (!any(test)) {
        stop("it seems that there are no events in the interval [", Tstart,
             ", ", Thoriz, ").\n")
    }
    id <- newdata[[id_var]]
    Time <- newdata[[Time_var]]
    event <- newdata[[event_var]]
    f <- factor(id, levels = unique(id))
    Time <- tapply(Time, f, tail, 1L)
    event <- tapply(event, f, tail, 1L)
    names(Time) <- names(event) <- as.character(unique(id))
    # subjects who had the event before Thoriz
    ind1 <- Time < Thoriz & event == 1
    # subjects who had the event after Thoriz
    ind2 <- Time > Thoriz
    # subjects who were censored in the interval (Tstart, Thoriz)
    ind3 <- Time < Thoriz & event == 0
    if (sum(ind1) < 5) {
        warning("there are fewer than 5 subjects with an event in the interval [",
                Tstart, ", ", Thoriz, ").\n")
    }
    tilde_Time <- pmin(Time, Thoriz)
    tilde_event <- event * as.numeric(Time < Thoriz)
    # newdata2 we set the event time at Tstart, and event to zero
    newdata2 <- newdata
    newdata2[[Time_var]] <- Tstart
    newdata2[[event_var]] <- 0
    # newdata2 we set the event time at tilde_Time
    newdata3 <- newdata2
    id. <- match(id, unique(id))
    ni <- tapply(id., id., length)
    newdata3[[Time_var]] <- rep(tilde_Time, ni)

    out <- if (!is_jm(object)) {
        # Super Learning
        V <- length(object) # number of folds
        L <- length(object[[1]]) # number of models
        tilde_Time_per_fold <-
            split(tilde_Time, tapply(newdata[["fold_"]], f, tail, 1L))
        run_over_folds <- function (v, object, newdata, newdata2, newdata3,
                                    tilde_Time, eps, L) {
            temp_q <- temp_q2 <- vector("list", L)
            for (l in seq_len(L)) {
                fold <- newdata2$fold_ == v
                # calculate Pr(T_i^* > \tilde T_i | T_i^* > t)
                preds <- predict(object[[v]][[l]], process = "event",
                                 times = tilde_Time[[v]],
                                 times_per_id = TRUE,
                                 newdata = newdata2[fold, ])
                pi_u_t <- preds$pred
                names(pi_u_t) <- preds$id
                # cumulative risk at tilde_Time
                f <- match(preds$id, unique(preds$id))
                pi_u_t <- tapply(pi_u_t, f, tail, n = 1L)
                # conditional survival probabilities
                temp_q[[l]] <- 1 - pi_u_t
                # calculate Pr(T_i^* > \tilde T_i + eps | T_i^* > \tilde T_i)
                preds2 <- predict(object[[v]][[l]], process = "event",
                                  times = tilde_Time[[v]] + eps,
                                  times_per_id = TRUE,
                                  newdata = newdata3[fold, ])
                pi_u_t2 <- preds2$pred
                names(pi_u_t2) <- preds2$id
                # cumulative risk at tilde_Time + eps
                f <- match(preds2$id, unique(preds2$id))
                pi_u_t2 <- tapply(pi_u_t2, f, tail, n = 1L)
                # conditional survival probabilities
                temp_q2[[l]] <- 1 - pi_u_t2
            }
            list(predictions = do.call("cbind", temp_q),
                 predictions2 = do.call("cbind", temp_q2))
        }
        cores <- min(cores, V)
        cl <- parallel::makeCluster(cores)
        invisible(parallel::clusterEvalQ(cl, library("JMbayes2")))
        res <-
            parallel::parLapply(cl, seq_len(V), run_over_folds, object = object,
                                newdata = newdata, newdata2 = newdata2,
                                newdata3 = newdata3,
                                tilde_Time = tilde_Time_per_fold, eps = eps,
                                L = L)
        parallel::stopCluster(cl)
        predictions <- do.call("rbind", lapply(res, "[[", "predictions"))
        predictions2 <- do.call("rbind", lapply(res, "[[", "predictions2"))
        weights_fun <- function (coefs, predictions, predictions2,
                                 tilde_event, eps) {
            coefs <- c(0.0, coefs)
            varpi <- exp(coefs) / sum(exp(coefs))
            qi_u_t <-
                rowSums(predictions * rep(varpi, each = nrow(predictions)))
            qi_u_t2 <-
                rowSums(predictions2 * rep(varpi, each = nrow(predictions2)))
            epce_fun(qi_u_t, qi_u_t2, tilde_event, eps)
        }
        opt <- optim(rep(0, L - 1), weights_fun, method = "BFGS",
                     predictions = predictions, predictions2 = predictions2,
                     tilde_event = tilde_event, eps = eps)
        coefs <- c(0, opt$par)
        varpi <- exp(coefs) / sum(exp(coefs))
        EPCE <- numeric(L)
        for (l in seq_len(L)) {
            EPCE[l] <- epce_fun(predictions[, l], predictions2[, l],
                                tilde_event, eps)
        }
        list(EPCE_per_model = EPCE, EPCE = opt$value, weights = varpi)
    } else {
        # calculate Pr(T_i^* > \tilde T_i | T_i^* > t)
        preds <- predict(object, newdata = newdata2, process = "event",
                         times = tilde_Time, times_per_id = TRUE)
        pi_u_t <- preds$pred
        names(pi_u_t) <- preds$id
        # cumulative risk at tilde_Time
        f <- match(preds$id, unique(preds$id))
        pi_u_t <- tapply(pi_u_t, f, tail, n = 1L)
        # conditional survival probabilities
        qi_u_t <- 1 - pi_u_t

        # calculate Pr(T_i^* > \tilde T_i + eps | T_i^* > \tilde T_i)
        preds2 <- predict(object, newdata = newdata3, process = "event",
                          times = tilde_Time + eps, times_per_id = TRUE)
        pi_u_t2 <- preds2$pred
        names(pi_u_t2) <- preds2$id
        # cumulative risk at tilde_Time
        f <- match(preds2$id, unique(preds2$id))
        pi_u_t2 <- tapply(pi_u_t2, f, tail, n = 1L)
        # conditional survival probabilities
        qi_u_t2 <- 1 - pi_u_t2
        # Calculate EPCE
        list(EPCE = epce_fun(qi_u_t, qi_u_t2, tilde_event, eps))
    }
    out$nr <- length(Time)
    out$nint <- sum(ind1)
    out$ncens <- sum(out$ind3)
    out$Tstart <- Tstart
    out$Thoriz <- Thoriz
    out$nfolds <- if (!is_jm(object)) length(object)
    out$nameObject <- deparse(substitute(object))
    class(out) <- "tvEPCE"
    out
}

print.tvEPCE <- function (x, digits = 4, ...) {
    if (!inherits(x, "tvEPCE"))
        stop("Use only with 'tvEPCE' objects.\n")
    if (!is.null(x$EPCE_per_model)) {
        cat("\nCross-Validated Expected Predictive Cross-Entropy using the Library of Joint Models '", x$nameObject, "'",
            sep = "")
            cat("\n\nSuper Learning Estimated EPCE:", round(x$EPCE, digits))
    } else {
        cat("\nExpected Predictive Cross-Entropy for the Joint Model '", x$nameObject, "'",
            sep = "")
        cat("\n\nEstimated EPCE:", round(x$EPCE, digits))
    }
    cat("\nIn the time interval: [", round(x$Tstart, digits),
            ", ", round(x$Thoriz, digits), ")", sep = "")
    cat("\nFor the ",  x$nr, " subjects at risk at time ",
        round(x$Tstart, digits), sep = "")
    cat("\nNumber of subjects with an event in [", round(x$Tstart, digits),
        ", ", round(x$Thoriz, digits), "): ", x$nint, sep = "")
    cat("\nNumber of subjects with a censored time in [", round(x$Tstart, digits),
        ", ", round(x$Thoriz, digits), "): ", x$ncens, sep = "")
    if (!is.null(x$EPCE_per_model)) {
        cat("\n\nBrier score per model:", round(x$EPCE_per_model, digits))
        cat("\nWeights per model:", round(x$weights, digits))
        cat("\nNumber of folds:", x$nfolds)
    }
    cat("\n\n")
    invisible(x)
}


tvBrier <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
                     integrated = FALSE, type_weights = c("model-based", "IPCW"),
                     cores = max(parallel::detectCores() - 1, 1), ...) {
    is_jm <- function (object) inherits(object, "jm")
    if (!is_jm(object)) {
        if (!all(sapply(unlist(object, recursive = FALSE), is_jm)))
            stop("Use only with 'jm' objects.\n")
    }
    if (is.null(Thoriz) && is.null(Dt)) {
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    }
    if (!is.null(Thoriz) && Thoriz <= Tstart) {
        stop("'Thoriz' must be larger than 'Tstart'.")
    }
    if (is.null(Thoriz)) {
        Thoriz <- Tstart + Dt
    }
    Tstart <- Tstart + 1e-06
    type_weights <- match.arg(type_weights)
    brier_fun <- function (pi_u_t, type_weights, weights,
                           ind1, ind2, ind3) {
        loss <- function (x) x * x
        res <- if (type_weights == "model-based") {
            events <- sum(loss(1 - pi_u_t[ind1]), na.rm = TRUE)
            no_events <- sum(loss(pi_u_t[ind2]), na.rm = TRUE)
            censored <- if (any(ind3)) {
                sum(weights * loss(1 - pi_u_t[ind3]) +
                        (1 - weights) * loss(pi_u_t[ind3]), na.rm = TRUE)
            } else 0.0
            (events + no_events + censored) / length(ind1)
        } else {
            mean(loss(as.numeric(ind1) - pi_u_t) * weights)
        }
        res
    }
    if (!is.data.frame(newdata) && !is.list(newdata) &&
        !all(sapply(newdata, is.data.frame))) {
        stop("'newdata' must be a data.frame or a list of data.frames.\n")
    }
    # if newdata is a list and not a data.frame,
    # Super Learning will be used
    if (!is.data.frame(newdata) && is.list(newdata)) {
        folds <- rep(seq_along(newdata), sapply(newdata, nrow))
        newdata <- do.call("rbind", newdata)
        newdata[["fold_"]] <- folds
    }
    # if Super Learning, object needs to be a list with length the
    # number of folds. In each element of the list, we have a list of fitted
    # models
    obj <- if (is_jm(object)) object else object[[1L]][[1L]]
    id_var <- obj$model_info$var_names$idVar
    time_var <- obj$model_info$var_names$time_var
    Time_var <- obj$model_info$var_names$Time_var
    event_var <- obj$model_info$var_names$event_var
    type_censoring <- object$model_info$type_censoring
    if (obj$model_info$CR_MS) {
        stop("'tvBrier()' currently only works for right censored data.")
    }
    if (is.null(newdata[[id_var]])) {
        stop("cannot find the '", id_var, "' variable in newdata.", sep = "")
    }
    if (is.null(newdata[[time_var]])) {
        stop("cannot find the '", time_var, "' variable in newdata.", sep = "")
    }
    if (any(sapply(Time_var, function (nmn) is.null(newdata[[nmn]])))) {
        stop("cannot find the '", paste(Time_var, collapse = ", "),
             "' variable(s) in newdata.", sep = "")
    }
    if (is.null(newdata[[event_var]])) {
        stop("cannot find the '", event_var, "' variable in newdata.", sep = "")
    }
    newdata <- newdata[newdata[[Time_var]] > Tstart, ]
    newdata <- newdata[newdata[[time_var]] <= Tstart, ]
    if (!nrow(newdata)) {
        stop("there are no data on subjects who had an observed event time after Tstart ",
             "and longitudinal measurements before Tstart.")
    }
    newdata[[id_var]] <- newdata[[id_var]][, drop = TRUE]

    br <- function (Thoriz) {
        test <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
        if (!any(test)) {
            stop("it seems that there are no events in the interval [", Tstart,
                 ", ", Thoriz, ").\n")
        }
        newdata2 <- newdata
        newdata2[[Time_var]] <- Tstart
        newdata2[[event_var]] <- 0

        id <- newdata[[id_var]]
        Time <- newdata[[Time_var]]
        event <- newdata[[event_var]]
        f <- factor(id, levels = unique(id))
        Time <- tapply(Time, f, tail, 1L)
        event <- tapply(event, f, tail, 1L)
        names(Time) <- names(event) <- as.character(unique(id))

        # subjects who had the event before Thoriz
        ind1 <- Time < Thoriz & event == 1
        # subjects who had the event after Thoriz
        ind2 <- Time > Thoriz
        # subjects who were censored in the interval (Tstart, Thoriz)
        ind3 <- Time < Thoriz & event == 0
        if (sum(ind1) < 5) {
            warning("there are fewer than 5 subjects with an event in the interval [",
                    Tstart, ", ", Thoriz, ").\n")
        }
        if (type_weights == "IPCW") {
            cens_data <- data.frame(Time = Time, cens_ind = 1 - event)
            censoring_dist <- survfit(Surv(Time, cens_ind) ~ 1, data = cens_data)
            weights <- numeric(length(Time))
            weights[ind1] <- 1 / summary(censoring_dist, times = Time[ind1])$surv
            weights[ind2] <- 1 / summary(censoring_dist, times = Thoriz)$surv
        }
        if (!is_jm(object)) {
            # Super Learning
            V <- length(object) # number of folds
            L <- length(object[[1]]) # number of models
            ids <- tapply(newdata2[[id_var]], newdata2[["fold_"]], unique)
            run_over_folds <- function (v, object, newdata, newdata2, type_weights,
                                        Tstart, Thoriz, ind1, ind2, ind3, ids, id,
                                        L) {
                temp_p <- temp_w <- vector("list", L)
                for (l in seq_len(L)) {
                    preds <- predict(object[[v]][[l]], process = "event",
                                     times = Thoriz,
                                     newdata = newdata2[newdata2$fold_ == v, ])
                    temp_p[[l]] <- preds$pred[preds$times > Tstart]
                    # which subjects in fold v had Time < Thoriz & event == 0
                    id_cens <- names(ind3[ind3])[names(ind3[ind3]) %in% ids[[v]]]
                    if (type_weights == "model-based" && length(id_cens)) {
                        preds2 <- predict(object[[v]][[l]],
                                          newdata = newdata[id %in% id_cens, ],
                                          process = "event", times = Thoriz)
                        weights <- preds2$pred
                        f <- factor(preds2$id, levels = unique(preds2$id))
                        names(weights) <- f
                        temp_w[[l]] <- tapply(weights, f, tail, 1)
                    }
                }
                list(predictions = do.call("cbind", temp_p),
                     W = if (type_weights == "model-based" && length(id_cens))
                         do.call("cbind", temp_w))
            }
            cores <- min(cores, V)
            cl <- parallel::makeCluster(cores)
            invisible(parallel::clusterEvalQ(cl, library("JMbayes2")))
            res <-
                parallel::parLapply(cl, seq_len(V), run_over_folds, object = object,
                                    newdata = newdata, newdata2 = newdata2,
                                    type_weights = type_weights, Tstart = Tstart,
                                    Thoriz = Thoriz, ind1 = ind1, ind2 = ind2,
                                    ind3 = ind3, ids = ids, id = id, L = L)
            parallel::stopCluster(cl)
            predictions <- do.call("rbind", lapply(res, "[[", "predictions"))
            W <- do.call("rbind", lapply(res, "[[", "W"))
            if (is.null(W)) {
                # two options: (i) IPCW, then W matrix of the weights
                # (ii) no censored observations, then matrix of zeros
                W <- matrix(if (type_weights == "IPCW") weights else 0.0,
                            length(weights), L)
            }
            list(predictions = predictions, W = W, ind1 = ind1, ind2 = ind2,
                 ind3 = ind3, Time = Time)
        } else {
            preds <- predict(object, newdata = newdata2, process = "event",
                             times = Thoriz)
            pi_u_t <- preds$pred
            names(pi_u_t) <- preds$id
            # cumulative risk at Thoriz
            pi_u_t <- pi_u_t[preds$times > Tstart]
            if (type_weights == "model-based" && any(ind3)) {
                nams <- names(ind3[ind3])
                preds2 <- predict(object, newdata = newdata[id %in% nams, ],
                                  process = "event", times = Thoriz)
                weights <- preds2$pred
                f <- factor(preds2$id, levels = unique(preds2$id))
                names(weights) <- f
                weights <- tapply(weights, f, tail, 1)
            }
            list(Brier = brier_fun(pi_u_t, type_weights, weights,
                                   ind1, ind2, ind3),
                 ind1 = ind1, ind2 = ind2, ind3 = ind3, Time = Time)
        }
    }
    out <- if (is_jm(object)) {
        if (integrated) {
            br1 <- br(0.5 * (Tstart + Thoriz))
            res <- br2 <- br(Thoriz)
            res$Brier <- 2 * br1$Brier / 3 + br2$Brier / 6
            res
        } else {
            br(Thoriz)
        }
    } else {
        temp <- if (integrated) {
            list(mid = br(0.5 * (Tstart + Thoriz)), last = br(Thoriz))
        } else {
            br(Thoriz)
        }
        weights_fun <- function (coefs, integrated, type_weights) {
            coefs <- c(0.0, coefs)
            varpi <- exp(coefs) / sum(exp(coefs))
            if (integrated) {
                ntemp <- length(temp)
                res <- numeric(ntemp)
                for (j in seq_len(ntemp)) {
                    tt <- temp[[j]]
                    pi_u_t <- rowSums(tt$predictions *
                                          rep(varpi, each = nrow(tt$predictions)))
                    weights <- if (type_weights == "model-based") {
                        rowSums(tt$W * rep(varpi, each = nrow(tt$W)))
                    } else tt$W
                    res[j] <- brier_fun(pi_u_t, type_weights, weights,
                                        tt$ind1, tt$ind2, tt$ind3)
                }
                2 * res[1L] / 3 + res[2L] / 6
            } else {
                pi_u_t <- rowSums(temp$predictions *
                                      rep(varpi, each = nrow(temp$predictions)))
                weights <- if (type_weights == "model-based") {
                    rowSums(temp$W * rep(varpi, each = nrow(temp$W)))
                } else temp$W
                brier_fun(pi_u_t, type_weights, weights, temp$ind1,
                          temp$ind2, temp$ind3)
            }
        }
        L <- length(object[[1]])
        opt <- optim(rep(0, L - 1), weights_fun, method = "BFGS",
                     integrated = integrated, type_weights = type_weights)
        coefs <- c(0, opt$par)
        varpi <- exp(coefs) / sum(exp(coefs))
        Brier <- numeric(L)
        for (l in seq_len(L)) {
            Brier[l] <- if (integrated) {
                tt_mid <- temp$mid
                br_mid <- brier_fun(tt_mid$predictions[, l], type_weights,
                                    tt_mid$W[, l], tt_mid$ind1, tt_mid$ind2,
                                    tt_mid$ind3)
                tt_last <- temp$last
                br_last <- brier_fun(tt_last$predictions[, l], type_weights,
                                     tt_last$W[, l], tt_last$ind1, tt_last$ind2,
                                     tt_last$ind3)
                2 * br_mid / 3 + br_last / 6
            } else {
                brier_fun(temp$predictions[, l], type_weights,
                                  temp$W[, l], temp$ind1, temp$ind2, temp$ind3)
            }
        }
        list(Brier = Brier, opt_Brier = opt$value, weights = varpi,
             Time = if (integrated) temp[[2]]$Time else temp$Time,
             ind1 = if (integrated) temp[[2]]$ind1 else temp$ind1,
             ind2 = if (integrated) temp[[2]]$ind2 else temp$ind2,
             ind3 = if (integrated) temp[[2]]$ind3 else temp$ind3)

    }
    out <- list(Brier = if (is_jm(object)) out$Brier else out$opt_Brier,
                Brier_per_model = if (!is_jm(object)) out$Brier,
                weights = if (!is_jm(object)) out$weights,
                nr = length(out$Time), nint = sum(out$ind1),
                ncens = sum(out$ind3), Tstart = Tstart, Thoriz = Thoriz,
                nfolds = if (!is_jm(object)) length(object),
                integrated = integrated, type_weights = type_weights,
                nameObject = deparse(substitute(object)))
    class(out) <- "tvBrier"
    out
}

print.tvBrier <- function (x, digits = 4, ...) {
    if (!inherits(x, "tvBrier"))
        stop("Use only with 'tvBrier' objects.\n")
    if (!is.null(x$Brier_per_model)) {
        cat("\nCross-Validated Prediction Error using the Library of Joint Models '", x$nameObject, "'",
            sep = "")
        if (x$integrated) {
            cat("\n\nSuper Learning Estimated Integrated Brier score:", round(x$Brier, digits))
        } else {
            cat("\n\nSuper Learning Estimated Brier score:", round(x$Brier, digits))
        }
    } else {
        cat("\nPrediction Error for the Joint Model '", x$nameObject, "'",
            sep = "")
        if (x$integrated) {
            cat("\n\nEstimated Integrated Brier score:", round(x$Brier, digits))
        } else {
            cat("\n\nEstimated Brier score:", round(x$Brier, digits))
        }
    }
    if (x$integrated) {
        cat("\nIn the time interval: [", round(x$Tstart, digits),
            ", ", round(x$Thoriz, digits), ")", sep = "")
    } else {
        cat("\nAt time:", round(x$Thoriz, digits))
    }
    cat("\nFor the ",  x$nr, " subjects at risk at time ",
        round(x$Tstart, digits), sep = "")
    cat("\nNumber of subjects with an event in [", round(x$Tstart, digits),
        ", ", round(x$Thoriz, digits), "): ", x$nint, sep = "")
    cat("\nNumber of subjects with a censored time in [", round(x$Tstart, digits),
        ", ", round(x$Thoriz, digits), "): ", x$ncens, sep = "")
    cat("\nAccounting for censoring using ",
        if (x$type_weights == "IPCW") "inverse probability of censoring Kaplan-Meier weights"
            else "model-based weights", sep = "")
    if (!is.null(x$Brier_per_model)) {
        cat("\n\nBrier score per model:", round(x$Brier_per_model, digits))
        cat("\nWeights per model:", round(x$weights, digits))
        cat("\nNumber of folds:", x$nfolds)
    }
    cat("\n\n")
    invisible(x)
}


CVdats <- create_folds(prothro, id_var = "id")
fit_models <- function (data) {
    library("JMbayes2")
    lmeFit <- lme(pro ~ time * treat, data = data, random = ~ time | id,
                  control = lmeControl(opt = "optim"))
    data_id <- data[!duplicated(data$id), ]
    CoxFit <- coxph(Surv(Time, death) ~ treat, data = data_id)
    jmFit1 <- jm(CoxFit, lmeFit, time_var = "time")
    jmFit2 <- update(jmFit1, functional_forms = ~ value(pro):treat)
    jmFit3 <- update(jmFit1, functional_forms = ~ area(pro))
    jmFit4 <- update(jmFit1, functional_forms = ~ area(pro):treat)
    jmFit5 <- update(jmFit1, functional_forms = ~ value(pro) + area(pro))
    out <- list(M1 = jmFit1, M2 = jmFit2, M3 = jmFit3, M4 = jmFit4, M5 = jmFit5)
    class(out) <- "jmList"
    out
}

cl <- parallel::makeCluster(5L)
Models_folds <- parallel::parLapply(cl, CVdats$training, fit_models)
parallel::stopCluster(cl)


tstr <- 3
thor <- 5
xxx1 <- tvBrier(Models_folds, CVdats$testing, integrated = TRUE,
                Tstart = tstr, Thoriz = thor)
xxx2 <- tvBrier(Models_folds, CVdats$testing, integrated = TRUE,
                type_weights = "IPCW", Tstart = tstr, Thoriz = thor)
xxx3 <- tvEPCE(Models_folds, CVdats$testing, Tstart = tstr, Thoriz = thor)

Models <- fit_models(prothro)
xxx4 <- tvEPCE(Models, prothro, Tstart = tstr, Thoriz = thor,
               model_weights = xxx3$weights)
xxx5 <- lapply(Models, tvEPCE, newdata = prothro, Tstart = tstr, Thoriz = thor)



ttt1 <- tvBrier(Models_folds[[1]][[4]], prothro, integrated = TRUE,
                Tstart = tstr, Thoriz = thor)
ttt2 <- tvBrier(Models_folds[[1]][[4]], prothro, integrated = TRUE,
                type_weights = "IPCW", Tstart = tstr, Thoriz = thor)
ttt3 <- tvEPCE(Models_folds[[1]][[4]], prothro, Tstart = tstr, Thoriz = thor)


xxx1 <- tvBrier(Models_folds, CVdats$testing, integrated = TRUE,
                Tstart = tstr, Thoriz = thor)
xxx2 <- tvBrier(Models_folds, CVdats$testing, integrated = TRUE,
                type_weights = "IPCW", Tstart = tstr, Thoriz = thor)



tstr <- 9
thor <- 10
ttt1 <- tvEPCE(Models_folds[[1]][[4]], aids, Tstart = tstr, Thoriz = thor)

yyy1 <- tvEPCE(Models_folds, CVdats$testing, Tstart = tstr, Thoriz = thor)



Models <- fit_models(prothro)
ND <- prothro[prothro$Time > tstr & prothro$time <= tstr, ]
ND$id <- ND$id[, drop = TRUE]
ND$Time <- tstr
ND$death <- 0
model_weights <- xxx1$weights

test <- predict(Models, model_weights, newdata = ND[ND$patient == 4, ],
                process = "event", return_newdata = TRUE)
test2 <- predict(Models, model_weights, newdata = ND[ND$patient == 4, ],
                 return_newdata = TRUE)
plot(test2, test)

test1 <- predict(Models[[5]], newdata = ND[ND$id == 3, ],
                 process = "event")

test2 <- predict(Models, weights = c(0, 0, 0, 0, 1), newdata = ND[ND$id == 3, ],
                 process = "event")

cbind(test1$pred, test2$pred)


#############################################################################
#############################################################################
#############################################################################

source("./Development/jm/PBC_data.R")
pbc2_lis <- split(pbc2, pbc2$id)
ind <- sample(length(pbc2_lis), 5)
pbc2_lis[ind] <- lapply(pbc2_lis[ind], function (d) {d$age <- NA; d})
ind <- sample(length(pbc2_lis), 5)
pbc2_lis[ind] <- lapply(pbc2_lis[ind], function (d) {d$sex <- NA; d})
pbc2 <- do.call("rbind", pbc2_lis)
pbc2.id <- pbc2[!duplicated(pbc2$id), ]
rm(ff, pbc2_lis, ind)

CVdats <- create_folds(pbc2, V = 5, id_var = "id")
fit_models <- function (data) {
    library("JMbayes2")
    data$status2 <- as.numeric(data$status != "alive")
    data_id <- data[!duplicated(data$id), ]
    lmeFit <- lme(log(serBilir) ~ year, data = data,
                  random = ~ year | id,
                  control = lmeControl(opt = "optim"))
    CoxFit <- coxph(Surv(years, status2) ~ 1, data = data_id)
    jmFit1 <- jm(CoxFit, lmeFit, time_var = "year")
    jmFit2 <- update(jmFit1,
                     functional_forms = ~ slope(log(serBilir)))
    jmFit3 <- update(jmFit1,
                     functional_forms = ~ value(log(serBilir)) + area(log(serBilir)))
    ###
    lmeFit2 <- lme(log(serBilir) ~ ns(year, 3, B = c(0, 14.4)) + sex + age,
                   data = data, random = ~ ns(year, 3, B = c(0, 14.4)) | id,
                   control = lmeControl(opt = "optim"))
    CoxFit2 <- coxph(Surv(years, status2) ~ sex + age, data = data_id)
    jmFit4 <- jm(CoxFit2, lmeFit2, time_var = "year")
    jmFit5 <- update(jmFit4,
                     functional_forms = ~ slope(log(serBilir)))
    out <- list(M1 = jmFit1, M2 = jmFit2, M3 = jmFit3, M4 = jmFit4, M5 = jmFit5)
    class(out) <- "jmList"
    out
}

cl <- parallel::makeCluster(5L)
Models_folds <- parallel::parLapply(cl, CVdats$training, fit_models)
parallel::stopCluster(cl)


tstr <- 6
thor <- 8

Brier_weights <- tvBrier(Models_folds, CVdats$testing, integrated = TRUE,
                         Tstart = tstr, Thoriz = thor)
Brier_weights

EPCE_weights <- tvEPCE(Models_folds, CVdats$testing, Tstart = tstr, Thoriz = thor)
EPCE_weights

Models <- fit_models(prothro)

ND <- pbc2[pbc2$years > tstr & pbc2$year <= tstr, ]
ND$id <- ND$id[, drop = TRUE]
ND$years <- tstr
ND$status2 <- 0


model_weights <- EPCE_weights$weights

predsEvent <- predict(Models, weights = model_weights, newdata = ND[ND$id == 8, ],
                      process = "event", return_newdata = TRUE)
predsLong <- predict(Models, weights = model_weights, newdata = ND[ND$id == 8, ],
                     return_newdata = TRUE)



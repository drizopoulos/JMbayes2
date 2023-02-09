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

CVdats <- create_folds(pbc2)
fit_models <- function (data) {
    library("JMbayes2")
    lmeFit <- lme(log(serBilir) ~ year * sex, data = data,
                  random = ~ year | id)
    data_id <- data[!duplicated(data$id), ]
    CoxFit <- coxph(Surv(years, status2) ~ sex, data = data_id)
    jmFit1 <- jm(CoxFit, lmeFit, time_var = "year")
    jmFit2 <- jm(CoxFit, lmeFit, time_var = "year",
                 functional_forms = ~ slope(log(serBilir)))
    jmFit3 <- jm(CoxFit, lmeFit, time_var = "year",
                 functional_forms = ~ area(log(serBilir)) + slope(log(serBilir)))
    list(M1 = jmFit1, M2 = jmFit2, M3 = jmFit3)
}

cl <- parallel::makeCluster(5L)
Models <- parallel::parLapply(cl, CVdats$training, fit_models)
parallel::stopCluster(cl)

if (FALSE) {
    object = Models[[1]][[1]]
    Tstart = 6
    Thoriz = 7
    cores = max(parallel::detectCores() - 1, 1)
    type_weights = "IPCW"
    newdata = pbc2

    object = Models
    Tstart = 6
    Thoriz = 7
    cores = max(parallel::detectCores() - 1, 1)
    type_weights = "IPCW"
    newdata = CVdats
}

tvBrier <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
                     type_weights = c("model-based", "IPCW"),
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
    Tstart <- Tstart
    Thoriz <- Thoriz
    type_weights <- match.arg(type_weights)
    brier_fun <- function (pi_u_t, type_weights, weights, ind1, ind2, ind3) {
        loss <- function (x) x * x
        if (type_weights == "model-based") {
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
    }
    if (!is.data.frame(newdata) &&
        (!is.list(newdata) && !all(names(newdata) %in% c("training", "testing")))) {
        stop("'newdata' must be a data.frame with more than one rows.\n")
    }
    # if newdata is a list with components 'training' and 'testing',
    # Super Learning will be used
    if (!is.data.frame(newdata) &&
        all(names(newdata) %in% c("training", "testing"))) {
        CV_data <- newdata
        newdata <- do.call("rbind", CV_data$testing)
        newdata[["fold_"]] <- rep(seq_along(CV_data$testing),
                                  sapply(CV_data$testing, nrow))
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
    test <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
    if (!any(test)) {
        stop("it seems that there are no events in the interval [Tstart, Thoriz).")
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
        warning("there are fewer than 5 subjects with an event in the interval ",
                "[Tstart, Thoriz).\n")
    }

    if (type_weights == "IPCW") {
        cens_data <- data.frame(Time = Time, cens_ind = 1 - event)
        censoring_dist <- survfit(Surv(Time, cens_ind) ~ 1, data = cens_data)
        weights <- numeric(length(Time))
        weights[ind1] <- 1 / summary(censoring_dist, times = Time[ind1])$surv
        weights[ind2] <- 1 / summary(censoring_dist, times = Thoriz)$surv
    }

    out <- if (!is_jm(object)) {
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
        if (is.null(W)) W <- matrix(weights, length(weights), L)
        weights_fun <- function (coefs, type_weights) {
            coefs <- c(0.0, coefs)
            varpi <- exp(coefs) / sum(exp(coefs))
            pi_u_t <- rowSums(predictions * rep(varpi, each = nrow(predictions)))
            names(pi_u_t) <- names(Time)
            if (type_weights == "model-based") {
                weights <- rowSums(W * rep(varpi, each = nrow(W)))
            }
            brier_fun(pi_u_t, type_weights, weights, ind1, ind2, ind3)
        }
        opt <- optim(rep(0, L - 1), weights_fun, method = "BFGS",
                     type_weights = type_weights)
        coefs <- c(0, opt$par)
        varpi <- exp(coefs) / sum(exp(coefs))
        Brier <- numeric(L)
        for (l in seq_len(L)) {
            Brier[l] <- brier_fun(predictions[, l], type_weights, W[, l],
                                  ind1, ind2, ind3)
        }
        list(Brier = Brier, opt_Brier = opt$value, weights = varpi)
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
        brier_fun(pi_u_t, type_weights, weights, ind1, ind2, ind3)
    }
    out <- list(Brier = if (is_jm(object)) out else out$opt_Brier,
                Brier_per_model = if (!is_jm(object)) out$Brier,
                weights = if (!is_jm(object)) out$weights,
                nr = length(Time), nint = sum(ind1),
                Tstart = Tstart, Thoriz = Thoriz,
                type_weights = type_weights,
                nameObject = deparse(substitute(object)))
    class(out) <- "tvBrier"
    out

}


print.tvBrier <- function (x, digits = 4, ...) {
    if (!inherits(x, "tvBrier"))
        stop("Use only with 'tvBrier' objects.\n")
    if (!is.null(x$Brier_per_model)) {
        cat("\nPrediction Error using the Joint Models '", x$nameObject, "'",
            sep = "")
        cat("\n\nSuper Learning Estimated Brier score:", round(x$Brier, digits))
    } else {
        cat("\nPrediction Error for the Joint Model '", x$nameObject, "'",
            sep = "")
        cat("\n\nEstimated Brier score:", round(x$Brier, digits))
    }
    cat("\nAt time:", round(x$Thoriz, digits))
    cat("\nUsing longitudinal information up to time: ", round(x$Tstart, digits),
        " (", x$nr, " subjects still at risk)", sep = "")
    cat("\nNumber of subjects with an event in the time interval [", round(x$Tstart, digits),
        ", ", round(x$Thoriz, digits), "): ", x$nint, sep = "")
    cat("\nAccounting for censoring using: ",
        if (x$type_weights == "IPCW") "inverse probability of censoring Kaplan-Meier weights"
            else "model-based weights", sep = "")
    if (!is.null(x$Brier_per_model)) {
        cat("\n\nBrier score per model:", round(x$Brier_per_model, digits))
        cat("\nweights per model:", round(x$weights, digits))
    }
    cat("\n\n")
    invisible(x)
}

ttt1 <- tvBrier(Models[[1]][[1]], pbc2, Tstart = 5, Thoriz = 7)
ttt2 <- tvBrier(Models[[1]][[1]], pbc2, type_weights = "IP", Tstart = 5, Thoriz = 7)

xxx1 <- tvBrier(Models, CVdats, Tstart = 5, Thoriz = 7)
xxx2 <- tvBrier(Models, CVdats, type_weights = "IP", Tstart = 5, Thoriz = 7)





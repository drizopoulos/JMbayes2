tvROC <- function (object, newdata, Tstart, ...) {
    UseMethod("tvROC")
}

tvROC.jm <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL, ...) {
    if (!inherits(object, "jm"))
        stop("Use only with 'jm' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (!is.null(Thoriz) && Thoriz <= Tstart)
        stop("'Thoriz' must be larger than 'Tstart'.")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    type_censoring <- object$model_info$type_censoring
    if (object$model_info$CR_MS)
        stop("'tvROC()' currently only works for right censored data.")
    Tstart <- Tstart + 1e-06
    Thoriz <- Thoriz + 1e-06
    id_var <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    Time_var <- object$model_info$var_names$Time_var
    event_var <- object$model_info$var_names$event_var
    if (is.null(newdata[[id_var]]))
        stop("cannot find the '", id_var, "' variable in newdata.", sep = "")
    if (is.null(newdata[[time_var]]))
        stop("cannot find the '", time_var, "' variable in newdata.", sep = "")
    if (any(sapply(Time_var, function (nmn) is.null(newdata[[nmn]]))))
        stop("cannot find the '", paste(Time_var, collapse = ", "),
             "' variable(s) in newdata.", sep = "")
    if (is.null(newdata[[event_var]]))
        stop("cannot find the '", event_var, "' variable in newdata.", sep = "")
    tt <- if (type_censoring == "right") newdata[[Time_var]] else newdata[[Time_var[2L]]]
    newdata[[id_var]] <- newdata[[id_var]][, drop = TRUE]
    id <- newdata[[id_var]]
    id <- match(id, unique(id))
    tt <- ave(tt, id, FUN = function (t) rep(tail(t, 1L) > Tstart, length(t)))
    newdata <- newdata[as.logical(tt), ]
    newdata <- newdata[newdata[[time_var]] <= Tstart, ]
    if (!nrow(newdata))
        stop("there are no data on subjects who had an observed event time after Tstart ",
             "and longitudinal measurements before Tstart.")
    test1 <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
    if (!any(test1))
        stop("it seems that there are no events in the interval [Tstart, Thoriz).")
    newdata2 <- newdata
    newdata2[[Time_var]] <- Tstart
    newdata2[[event_var]] <- 0
    preds <- predict(object, newdata = newdata2, process = "event",
                     times = Thoriz, ...)
    qi_u_t <- 1 - preds$pred
    names(qi_u_t) <- preds$id
    qi_u_t <- qi_u_t[preds$times > Tstart]

    id <- newdata[[id_var]]
    Time <- newdata[[Time_var]]
    event <- newdata[[event_var]]
    f <- factor(id, levels = unique(id))
    Time <- tapply(Time, f, tail, 1L)
    event <- tapply(event, f, tail, 1L)
    names(Time) <- names(event) <- as.character(unique(id))
    # subjects who died before Thoriz
    ind1 <- Time < Thoriz & event == 1
    # subjects who were censored in the interval (Tstart, Thoriz)
    ind2 <- Time < Thoriz & event == 0
    ind <- ind1 | ind2
    if (any(ind2)) {
        nams <- names(ind2[ind2])
        preds2 <- predict(object, newdata = newdata[id %in% nams, ],
                          process = "event", times = Thoriz, ...)
        pi_u_t <- preds2$pred
        f <- factor(preds2$id, levels = unique(preds2$id))
        names(pi_u_t) <- f
        pi_u_t <- tapply(pi_u_t, f, tail, 1)
        nams2 <- names(ind2[ind2])
        ind[ind2] <- ind[ind2] * pi_u_t[nams2]
    }
    # calculate sensitivity and specificity
    thrs <- seq(0, 1, length = 101)
    Check <- outer(qi_u_t, thrs, "<")
    nTP <- colSums(Check * c(ind))
    nFN <- sum(ind) - nTP
    TP <- nTP / sum(ind)
    nFP <- colSums(Check * c(1 - ind))
    nTN <- sum(1 - ind) - nFP
    FP <- nFP / sum(1 - ind)
    Q <- colMeans(Check)
    Q. <- 1 - Q
    k.1.0 <- (TP - Q) / Q.
    k.0.0 <- (1 - FP - Q.) / Q
    P <- mean(ind)
    P. <- 1 - P
    k.05.0 <- (P * Q. * k.1.0 + P. * Q * k.0.0) / (P * Q. + P. * Q)
    f1score <- 2 * nTP / (2 * nTP + nFN + nFP)
    F1score <- median(thrs[f1score == max(f1score)])
    youden <- TP - FP
    Youden <- median(thrs[youden == max(youden)])
    out <- list(TP = TP, FP = FP, nTP = nTP, nFN = nFN, nFP = nFP, nTN = nTN,
                qSN = k.1.0, qSP = k.0.0, qOverall = k.05.0,
                thrs = thrs, F1score = F1score, Youden = Youden,
                Tstart = Tstart, Thoriz = Thoriz, nr = length(unique(id)),
                classObject = class(object),
                nameObject = deparse(substitute(object)))
    class(out) <- "tvROC"
    out
}

print.tvROC <- function (x, digits = 4, ...) {
    cat("\n\tTime-dependent Sensitivity and Specificity for the Joint Model",
        x$nameObject)
    cat("\n\nAt time:", round(x$Thoriz, digits))
    cat("\nUsing information up to time: ", round(x$Tstart, digits),
        " (", x$nr, " subjects still at risk)\n\n", sep = "")
    d <- data.frame("cut-off" = x$thrs, "SN" = x$TP, "SP" = 1 - x$FP,
                    "qSN" = x$qSN, "qSP" = x$qSP, check.names = FALSE,
                    check.rows = FALSE)
    xx <- rep("", nrow(d))
    xx[which.min(abs(x$thr - x$Youden))] <- "*"
    d[[" "]] <- xx
    d <- d[!is.na(d$qSN) & !is.na(d$qSP), ]
    d <- d[!duplicated(d[c("SN", "SP")]), ]
    row.names(d) <- 1:nrow(d)
    print(d)
    cat("\n")
    invisible(x)
}

plot.tvROC <- function (x, legend = FALSE,
                        optimal_cutoff = c("", "F1score", "Youden"),
                        xlab = "1 - Specificity", ylab = "Sensitivity", ...) {
    plot(x$FP, x$TP, type = "l", xlab = xlab, ylab = ylab, ...)
    abline(a = 0, b = 1, lty = 3)
    optimal_cutoff <- match.arg(optimal_cutoff)
    if (optimal_cutoff == "F1score")
        abline(v = x$thrs[which.max(x$F1score)], lty = 3, lwd = 2, col = 2)
    if (optimal_cutoff == "Youden")
        abline(v = x$thrs[which.max(x$TP - x$FP)], lty = 3, lwd = 2, col = 2)
    if (legend) {
        legend("bottomright", c(paste("At time:", round(x$Thoriz, 1), "\n"),
                                paste("Using information up to time:",
                                      round(x$Tstart, 1))),
               bty = "n")
    }
    invisible()
}

tvAUC <- function (object, newdata, Tstart, ...) {
    UseMethod("tvAUC")
}

tvAUC.jm <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL, ...) {
    if (!inherits(object, "jm"))
        stop("Use only with 'jm' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (!is.null(Thoriz) && Thoriz <= Tstart)
        stop("'Thoriz' must be larger than 'Tstart'.")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    type_censoring <- object$model_info$type_censoring
    if (object$model_info$CR_MS)
        stop("'tvROC()' currently only works for right censored data.")
    Tstart <- Tstart + 1e-06
    Thoriz <- Thoriz + 1e-06
    id_var <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    Time_var <- object$model_info$var_names$Time_var
    event_var <- object$model_info$var_names$event_var
    if (is.null(newdata[[id_var]]))
        stop("cannot find the '", id_var, "' variable in newdata.", sep = "")
    if (is.null(newdata[[time_var]]))
        stop("cannot find the '", time_var, "' variable in newdata.", sep = "")
    if (any(sapply(Time_var, function (nmn) is.null(newdata[[nmn]]))))
        stop("cannot find the '", paste(Time_var, collapse = ", "),
             "' variable(s) in newdata.", sep = "")
    if (is.null(newdata[[event_var]]))
        stop("cannot find the '", event_var, "' variable in newdata.", sep = "")
    newdata <- newdata[order(newdata[[Time_var]]), ]
    newdata <- newdata[newdata[[Time_var]] > Tstart, ]
    newdata <- newdata[newdata[[time_var]] <= Tstart, ]
    newdata[[id_var]] <- newdata[[id_var]][, drop = TRUE]
    test1 <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
    if (!any(test1))
        stop("it seems that there are no events in the interval [Tstart, Thoriz).")
    newdata2 <- newdata
    newdata2[[Time_var]] <- Tstart
    newdata2[[event_var]] <- 0
    preds <- predict(object, newdata = newdata2, process = "event",
                     times = Thoriz, ...)
    si_u_t <- 1 - preds$pred
    names(si_u_t) <- preds$id
    si_u_t <- si_u_t[preds$times > Tstart]

    id <- newdata[[id_var]]
    Time <- newdata[[Time_var]]
    event <- newdata[[event_var]]
    f <- factor(id, levels = unique(id))
    Time <- tapply(Time, f, tail, 1L)
    event <- tapply(event, f, tail, 1L)
    names(Time) <- names(event) <- as.character(unique(id))
    if (any(dupl <- duplicated(Time))) {
        Time[dupl] <- Time[dupl] + runif(length(Time[dupl]), 1e-07, 1e-06)
    }
    if (!all(names(si_u_t) == names(Time)))
        stop("mismatch between 'Time' variable names and survival probabilities names.")

    auc <- if (length(Time) > 1L) {
        pairs <- combn(as.character(unique(id)), 2)
        Ti <- Time[pairs[1, ]]
        Tj <- Time[pairs[2, ]]
        di <- event[pairs[1, ]]
        dj <- event[pairs[2, ]]
        si_u_t_i <- si_u_t[pairs[1, ]]
        si_u_t_j <- si_u_t[pairs[2, ]]
        ind1 <- (Ti <= Thoriz & di == 1) & Tj > Thoriz
        ind2 <- (Ti <= Thoriz & di == 0) & Tj > Thoriz
        ind3 <- (Ti <= Thoriz & di == 1) & (Tj <= Thoriz & dj == 0)
        ind4 <- (Ti <= Thoriz & di == 0) & (Tj <= Thoriz & dj == 0)
        names(ind1) <- names(ind2) <- names(ind3) <- names(ind4) <-
            paste(names(Ti), names(Tj), sep = "_")
        ind <- ind1 | ind2 | ind3 | ind4
        if (any(ind2)) {
            nams <- strsplit(names(ind2[ind2]), "_")
            nams_i <- sapply(nams, "[", 1)
            unq_nams_i <- unique(nams_i)
            preds2 <- predict(object, newdata = newdata[id %in% unq_nams_i, ],
                              process = "event", times = Thoriz, ...)
            pi_u_t <- preds2$pred
            f <- factor(preds2$id, levels = unique(preds2$id))
            names(pi_u_t) <- f
            pi_u_t <- tapply(pi_u_t, f, tail, 1)
            ind[ind2] <- ind[ind2] * pi_u_t[nams_i]
        }
        if (any(ind3)) {
            nams <- strsplit(names(ind3[ind3]), "_")
            nams_j <- sapply(nams, "[", 2)
            unq_nams_j <- unique(nams_j)
            preds3 <- predict(object, newdata = newdata[id %in% unq_nams_j, ],
                              process = "event", times = Thoriz)
            qi_u_t <- preds3$pred
            f <- factor(preds3$id, levels = unique(preds3$id))
            names(qi_u_t) <- f
            qi_u_t <- 1 - tapply(qi_u_t, f, tail, 1)
            ind[ind3] <- ind[ind3] * qi_u_t[nams_j]
        }
        if (any(ind4)) {
            nams <- strsplit(names(ind4[ind4]), "_")
            nams_i <- sapply(nams, "[", 1)
            nams_j <- sapply(nams, "[", 2)
            unq_nams_i <- unique(nams_i)
            unq_nams_j <- unique(nams_j)
            preds4_i <- predict(object, newdata = newdata[id %in% unq_nams_i, ],
                                process = "event", times = Thoriz, ...)
            pi_u_t <- preds4_i$pred
            f <- factor(preds4_i$id, levels = unique(preds4_i$id))
            names(pi_u_t) <- f
            pi_u_t <- tapply(pi_u_t, f, tail, 1)

            preds4_j <- predict(object, newdata = newdata[id %in% unq_nams_j, ],
                                process = "event", times = Thoriz, ...)
            qi_u_t <- preds4_j$pred
            f <- factor(preds4_j$id, levels = unique(preds4_j$id))
            names(qi_u_t) <- f
            qi_u_t <- 1 - tapply(qi_u_t, f, tail, 1)
            ind[ind4] <- ind[ind4] * pi_u_t[nams_i] * qi_u_t[nams_j]
        }
        sum((si_u_t_i < si_u_t_j) * c(ind), na.rm = TRUE) / sum(ind, na.rm = TRUE)
    } else {
        NA
    }
    out <- list(auc = auc, Tstart = Tstart, Thoriz = Thoriz, nr = length(unique(id)),
                classObject = class(object), nameObject = deparse(substitute(object)))
    class(out) <- "tvAUC"
    out
}

tvAUC.jm <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL, ...) {
    roc <- tvROC(object, newdata, Tstart, Thoriz, Dt, ...)
    TP <- roc$TP
    FP <- roc$FP
    auc <- sum(0.5 * diff(FP) * (TP[-1L] + TP[-length(TP)]), na.rm = TRUE)
    out <- list(auc = auc, Tstart = Tstart, Thoriz = roc$Thoriz, nr = roc$nr,
                classObject = class(object), nameObject = deparse(substitute(object)))
    class(out) <- "tvAUC"
    out
}

tvAUC.tvROC <- function (object, ...) {
    TP <- object$TP
    FP <- object$FP
    auc <- sum(0.5 * diff(FP) * (TP[-1L] + TP[-length(TP)]), na.rm = TRUE)
    out <- list(auc = auc, Tstart = object$Tstart, Thoriz = object$Thoriz,
                nr = object$nr, classObject = object$classObject,
                nameObject = object$nameObject)
    class(out) <- "tvAUC"
    out
}

print.tvAUC <- function (x, digits = 4, ...) {
    if (!inherits(x, "tvAUC"))
        stop("Use only with 'tvAUC' objects.\n")
    if (x$class == "jm")
        cat("\n\tTime-dependent AUC for the Joint Model",  x$nameObject)
    else
        cat("\n\tTime-dependent AUC for the Cox Model",  x$nameObject)
    cat("\n\nEstimated AUC:", round(x$auc, digits))
    cat("\nAt time:", round(x$Thoriz, digits))
    cat("\nUsing information up to time: ", round(x$Tstart, digits),
        " (", x$nr, " subjects still at risk)", sep = "")
    cat("\n\n")
    invisible(x)
}

calibration_plot <- function (object, newdata, Tstart, Thoriz = NULL,
                              Dt = NULL, df_ns = 3, plot = TRUE,
                              add_density = TRUE, col = "red", lty = 1, lwd = 1,
                              col_dens = "grey",
                              xlab = "Predicted Probabilities",
                              ylab = "Observed Probabilities", main = "", ...) {
    if (!inherits(object, "jm"))
        stop("Use only with 'jm' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (!is.null(Thoriz) && Thoriz <= Tstart)
        stop("'Thoriz' must be larger than 'Tstart'.")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    type_censoring <- object$model_info$type_censoring
    if (object$model_info$CR_MS)
        stop("'tvROC()' currently only works for right censored data.")
    Tstart <- Tstart + 1e-06
    Thoriz <- Thoriz + 1e-06
    id_var <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    Time_var <- object$model_info$var_names$Time_var
    event_var <- object$model_info$var_names$event_var
    if (is.null(newdata[[id_var]]))
        stop("cannot find the '", id_var, "' variable in newdata.", sep = "")
    if (is.null(newdata[[time_var]]))
        stop("cannot find the '", time_var, "' variable in newdata.", sep = "")
    if (any(sapply(Time_var, function (nmn) is.null(newdata[[nmn]]))))
        stop("cannot find the '", paste(Time_var, collapse = ", "),
             "' variable(s) in newdata.", sep = "")
    if (is.null(newdata[[event_var]]))
        stop("cannot find the '", event_var, "' variable in newdata.", sep = "")
    newdata <- newdata[newdata[[Time_var]] > Tstart, ]
    newdata <- newdata[newdata[[time_var]] <= Tstart, ]
    if (!nrow(newdata))
        stop("there are no data on subjects who had an observed event time after Tstart ",
             "and longitudinal measurements before Tstart.")
    newdata[[id_var]] <- newdata[[id_var]][, drop = TRUE]
    test1 <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
    if (!any(test1))
        stop("it seems that there are no events in the interval [Tstart, Thoriz).")
    newdata2 <- newdata
    newdata2[[Time_var]] <- Tstart
    newdata2[[event_var]] <- 0
    preds <- predict(object, newdata = newdata2, process = "event",
                     times = Thoriz, ...)
    pi_u_t <- preds$pred
    names(pi_u_t) <- preds$id
    pi_u_t <- pi_u_t[preds$times > Tstart]

    id <- newdata[[id_var]]
    Time <- newdata[[Time_var]]
    event <- newdata[[event_var]]
    f <- factor(id, levels = unique(id))
    Time <- tapply(Time, f, tail, 1L)
    event <- tapply(event, f, tail, 1L)
    names(Time) <- names(event) <- as.character(unique(id))
    cal_DF <- data.frame(Time = Time, event = event, preds = pi_u_t[names(Time)])
    cloglog <- function (x) log(-log(1 - x))
    Bounds <- quantile(cloglog(pi_u_t), probs = c(0.05, 0.95))
    form <- paste0("ns(cloglog(preds), df = ", df_ns,
                   ", B = c(", round(Bounds[1L], 2), ", ",
                   round(Bounds[2L], 2), "))")
    form <- paste("Surv(Time, event) ~", form)
    cal_Cox <- coxph(as.formula(form), data = cal_DF)
    qs <- quantile(pi_u_t, probs = c(0.01, 0.99))
    probs_grid <- data.frame(preds = seq(qs[1L], qs[2L], length.out = 100L))
    obs <- 1 - c(summary(survfit(cal_Cox, newdata = probs_grid), times = Thoriz)$surv)
    obs_pi_u_t <- 1 - c(summary(survfit(cal_Cox, newdata = cal_DF), times = Thoriz)$surv)
    if (plot) {
        plot(probs_grid$preds, obs, type = "l", col = col, lwd = lwd, lty = lty,
             xlab = xlab, ylab = ylab, main = main, xlim = c(0, 1),
             ylim = c(0, 1))
        abline(0, 1, lty = 2)
        if (add_density) {
            par(new = TRUE)
            plot(density(pi_u_t), axes = FALSE, xlab = "", ylab = "", main = "",
                 col = col_dens)
            axis(side = 4)
        }
        invisible()
    } else {
        list("observed" = obs, "predicted" = probs_grid$preds,
             "pi_u_t" = pi_u_t, "obs_pi_u_t" = obs_pi_u_t)
    }
}

calibration_metrics <- function (object, newdata, Tstart, Thoriz = NULL,
                                 Dt = NULL, df_ns = 3, ...) {
    comps <- calibration_plot(object, newdata, Tstart, Thoriz = Thoriz, Dt = Dt,
                              df_ns = df_ns, plot = FALSE)
    diff <- abs(as.vector(comps$pi_u_t) - comps$obs_pi_u_t)
    ICI <- mean(diff)
    E50 <- median(diff)
    E90 <- unname(quantile(diff, probs = 0.9))
    c("ICI" = ICI, "E50" = E50, "E90" = E90)
}

tvEPCE <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
                    eps = 0.001, model_weights = NULL,
                    parallel = c("snow", "multicore"),
                    cores = parallelly::availableCores(omit = 1L), ...) {
    parallel <- match.arg(parallel)
    is_jm <- function (object) inherits(object, "jm")
    is_jmList <- function (object) inherits(object, "jmList")
    if (!is_jm(object)) {
        if (!all(sapply(unlist(object, recursive = FALSE), is_jm)) &&
            !all(sapply(object, is_jm)))
            stop("Use only with 'jm' objects.\n")
        if (!is.null(model_weights) && !is_jmList(object)) {
            stop("When 'model_weights' is not NULL, 'object' must have the class ",
                 "'jmList'.\n")
        }
        if (is_jmList(object) && is.null(model_weights)) {
            stop("For 'jmList' objects 'model_weights' cannot be NULL.\n")
        }
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
        SL <- TRUE
        folds <- rep(seq_along(newdata), sapply(newdata, nrow))
        newdata <- do.call("rbind", newdata)
        newdata[["fold_"]] <- folds
    } else SL <- FALSE
    # if Super Learning, object needs to be a list with length the
    # number of folds. In each element of the list, we have a list of fitted
    # models
    obj <- if (is_jm(object)) object else if (is_jmList(object)) object[[1L]] else object[[1L]][[1L]]
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
    out <- if (!is_jm(object) && SL) {
        # Super Learning
        V <- length(object) # number of folds
        L <- length(object[[1]]) # number of models
        tilde_Time_per_fold <-
            split(tilde_Time, tapply(newdata[["fold_"]], f, tail, 1L))
        run_over_folds <- function (v, object, newdata, newdata2, newdata3,
                                    tilde_Time, eps, L, parallel, cores = 1L) {
            temp_q <- temp_q2 <- vector("list", L)
            for (l in seq_len(L)) {
                fold <- newdata2$fold_ == v
                # calculate Pr(T_i^* > \tilde T_i | T_i^* > t)
                preds <- predict(object[[v]][[l]], process = "event",
                                 times = tilde_Time[[v]],
                                 times_per_id = TRUE, parallel = parallel,
                                 cores = cores, newdata = newdata2[fold, ])
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
                                  times_per_id = TRUE, parallel = parallel,
                                  cores = cores, newdata = newdata3[fold, ])
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
        if (cores > 1L) {
            have_mc <- have_snow <- FALSE
            if (parallel == "multicore") {
                have_mc <- .Platform$OS.type != "windows"
            } else if (parallel == "snow") {
                have_snow <- TRUE
            }
            if (!have_mc && !have_snow) cores <- 1L
            loadNamespace("parallel")
        }
        if (cores > 1L) {
            cores2 <- 1 # parallelly::availableCores(omit = 1)
            if (have_mc) {
                res <-
                    parallel::mclapply(seq_len(V), run_over_folds, object = object,
                                       newdata = newdata, newdata2 = newdata2,
                                       newdata3 = newdata3,
                                       tilde_Time = tilde_Time_per_fold, eps = eps,
                                       L = L, parallel = parallel, cores = cores2,
                                       mc.cores = cores)
            } else {
                cl <- parallel::makePSOCKcluster(rep("localhost", cores))
                invisible(parallel::clusterEvalQ(cl, library("JMbayes2")))
                res <-
                    parallel::parLapply(cl, seq_len(V), run_over_folds, object = object,
                                        newdata = newdata, newdata2 = newdata2,
                                        newdata3 = newdata3,
                                        tilde_Time = tilde_Time_per_fold, eps = eps,
                                        L = L, parallel = parallel, cores = cores2)
                parallel::stopCluster(cl)
            }
        } else {
            res <- lapply(seq_len(V), run_over_folds, object = object,
                          newdata = newdata, newdata2 = newdata2,
                          newdata3 = newdata3,
                          tilde_Time = tilde_Time_per_fold, eps = eps,
                          L = L, parallel = parallel)
        }
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
        preds <- if (is_jm(object)) {
            predict(object, newdata = newdata2, process = "event",
                    times = tilde_Time, times_per_id = TRUE,
                    parallel = parallel)
        } else if (is_jmList(object)) {
            predict(object, newdata = newdata2, process = "event",
                    times = tilde_Time, times_per_id = TRUE,
                    parallel = parallel, weights = model_weights)
        }
        pi_u_t <- preds$pred
        names(pi_u_t) <- preds$id
        # cumulative risk at tilde_Time
        f <- match(preds$id, unique(preds$id))
        pi_u_t <- tapply(pi_u_t, f, tail, n = 1L)
        # conditional survival probabilities
        qi_u_t <- 1 - pi_u_t

        # calculate Pr(T_i^* > \tilde T_i + eps | T_i^* > \tilde T_i)
        preds2 <- if (is_jm(object)) {
            predict(object, newdata = newdata3, process = "event",
                    times = tilde_Time + eps, times_per_id = TRUE,
                    parallel = parallel)
        } else if (is_jmList(object)) {
            predict(object, newdata = newdata3, process = "event",
                    times = tilde_Time + eps, times_per_id = TRUE,
                    parallel = parallel, weights = model_weights)
        }
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
    out$ncens <- sum(ind3)
    out$Tstart <- Tstart
    out$Thoriz <- Thoriz
    out$nfolds <- if (!is_jm(object) && !is_jmList(object)) length(object)
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
        cat("\n\nEPCE per model:", round(x$EPCE_per_model, digits))
        cat("\nWeights per model:", round(x$weights, digits))
        cat("\nNumber of folds:", x$nfolds)
    }
    cat("\n\n")
    invisible(x)
}


tvBrier <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
                     integrated = FALSE, type_weights = c("model-based", "IPCW"),
                     model_weights = NULL, parallel = c("snow", "multicore"),
                     cores = parallelly::availableCores(omit = 1L), ...) {
    parallel <- match.arg(parallel)
    is_jm <- function (object) inherits(object, "jm")
    is_jmList <- function (object) inherits(object, "jmList")
    if (!is_jm(object)) {
        if (!all(sapply(unlist(object, recursive = FALSE), is_jm)) &&
            !all(sapply(object, is_jm)))
            stop("Use only with 'jm' objects.\n")
        if (!is.null(model_weights) && !is_jmList(object)) {
            stop("When 'model_weights' is not NULL, 'object' must have the class ",
                 "'jmList'.\n")
        }
        if (is_jmList(object) && is.null(model_weights)) {
            stop("For 'jmList' objects 'model_weights' cannot be NULL.\n")
        }
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
        SL <- TRUE
        folds <- rep(seq_along(newdata), sapply(newdata, nrow))
        newdata <- do.call("rbind", newdata)
        newdata[["fold_"]] <- folds
    } else SL <- FALSE
    # if Super Learning, object needs to be a list with length the
    # number of folds. In each element of the list, we have a list of fitted
    # models
    obj <- if (is_jm(object)) object else if (is_jmList(object)) object[[1L]] else object[[1L]][[1L]]
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
        if (!is_jm(object) && SL) {
            # Super Learning
            V <- length(object) # number of folds
            L <- length(object[[1]]) # number of models
            ids <- tapply(newdata2[[id_var]], newdata2[["fold_"]], unique)
            run_over_folds <- function (v, object, newdata, newdata2, type_weights,
                                        Tstart, Thoriz, ind1, ind2, ind3, ids, id,
                                        L, parallel, cores = 1L) {
                temp_p <- temp_w <- vector("list", L)
                for (l in seq_len(L)) {
                    preds <- predict(object[[v]][[l]], process = "event",
                                     times = Thoriz, parallel = parallel,
                                     cores = cores,
                                     newdata = newdata2[newdata2$fold_ == v, ])
                    temp_p[[l]] <- preds$pred[preds$times > Tstart]
                    # which subjects in fold v had Time < Thoriz & event == 0
                    id_cens <- names(ind3[ind3])[names(ind3[ind3]) %in% ids[[v]]]
                    if (type_weights == "model-based" && length(id_cens)) {
                        preds2 <- predict(object[[v]][[l]],
                                          newdata = newdata[id %in% id_cens, ],
                                          process = "event", times = Thoriz,
                                          parallel = parallel, cores = cores)
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
            if (cores > 1L) {
                have_mc <- have_snow <- FALSE
                if (parallel == "multicore") {
                    have_mc <- .Platform$OS.type != "windows"
                } else if (parallel == "snow") {
                    have_snow <- TRUE
                }
                if (!have_mc && !have_snow) cores <- 1L
                loadNamespace("parallel")
            }
            if (cores > 1L) {
                cores2 <- 1 #parallelly::availableCores()
                if (have_mc) {
                    res <-
                        parallel::mclapply(seq_len(V), run_over_folds, object = object,
                                           newdata = newdata, newdata2 = newdata2,
                                           type_weights = type_weights, Tstart = Tstart,
                                           Thoriz = Thoriz, ind1 = ind1, ind2 = ind2,
                                           ind3 = ind3, ids = ids, id = id, L = L,
                                           parallel = parallel, cores = cores2,
                                           mc.cores = cores)
                } else {
                    cl <- parallel::makePSOCKcluster(rep("localhost", cores))
                    invisible(parallel::clusterEvalQ(cl, library("JMbayes2")))
                    res <-
                        parallel::parLapply(cl, seq_len(V), run_over_folds, object = object,
                                            newdata = newdata, newdata2 = newdata2,
                                            type_weights = type_weights, Tstart = Tstart,
                                            Thoriz = Thoriz, ind1 = ind1, ind2 = ind2,
                                            ind3 = ind3, ids = ids, id = id, L = L,
                                            parallel = parallel, cores = cores2)
                    parallel::stopCluster(cl)
                }
            } else {
                res <-
                    lapply(seq_len(V), run_over_folds, object = object,
                           newdata = newdata, newdata2 = newdata2,
                           type_weights = type_weights, Tstart = Tstart,
                           Thoriz = Thoriz, ind1 = ind1, ind2 = ind2,
                           ind3 = ind3, ids = ids, id = id, L = L,
                           parallel = parallel)
            }
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
            preds <- if (is_jm(object)) {
                predict(object, newdata = newdata2, process = "event",
                        times = Thoriz, parallel = parallel)
            } else if (is_jmList(object)) {
                predict(object, newdata = newdata2, process = "event",
                        times = Thoriz, parallel = parallel,
                        weights = model_weights)
            }
            pi_u_t <- preds$pred
            names(pi_u_t) <- preds$id
            # cumulative risk at Thoriz
            pi_u_t <- pi_u_t[preds$times > Tstart]
            if (type_weights == "model-based" && any(ind3)) {
                nams <- names(ind3[ind3])
                preds2 <- if (is_jm(object)) {
                    predict(object, newdata = newdata[id %in% nams, ],
                            process = "event", times = Thoriz,
                            parallel = parallel)
                } else if (is_jmList(object)) {
                    predict(object, newdata = newdata[id %in% nams, ],
                            process = "event", times = Thoriz,
                            parallel = parallel, weights = model_weights)
                }
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
    out <- if (is_jm(object) || is_jmList(object)) {
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
    out <- list(Brier = if (is_jm(object) || is_jmList(object)) out$Brier else out$opt_Brier,
                Brier_per_model = if (!is_jm(object)) out$Brier,
                weights = if (!is_jm(object) && !is_jmList(object)) out$weights,
                nr = length(out$Time), nint = sum(out$ind1),
                ncens = sum(out$ind3), Tstart = Tstart, Thoriz = Thoriz,
                nfolds = if (!is_jm(object) && !is_jmList(object)) length(object),
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
        if (x$integrated) {
            cat("\n\nIntegrated Brier score per model:", round(x$Brier_per_model, digits))
        } else {
            cat("\n\nBrier score per model:", round(x$Brier_per_model, digits))
        }
        cat("\nWeights per model:", round(x$weights, digits))
        cat("\nNumber of folds:", x$nfolds)
    }
    cat("\n\n")
    invisible(x)
}

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

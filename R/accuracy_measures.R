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
        stop("there are no data on subjects who had an observed event time after Tstart",
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
        stop("there are no data on subjects who had an observed event time after Tstart",
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
    comps <- calibration_plot(object, newdata, Tstart, Dt = Dt, df_ns = df_ns,
                              plot = FALSE)
    diff <- abs(as.vector(comps$pi_u_t) - comps$obs_pi_u_t)
    ICI <- mean(diff)
    E50 <- median(diff)
    E90 <- unname(quantile(diff, probs = 0.9))
    c("ICI" = ICI, "E50" = E50, "E90" = E90)
}

tvBrier <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL, ...) {
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
        stop("there are no data on subjects who had an observed event time after Tstart",
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
    # cumulative risk at Thoriz
    pi_u_t <- pi_u_t[preds$times > Tstart]

    id <- newdata[[id_var]]
    Time <- newdata[[Time_var]]
    event <- newdata[[event_var]]
    f <- factor(id, levels = unique(id))
    Time <- tapply(Time, f, tail, 1L)
    event <- tapply(event, f, tail, 1L)
    names(Time) <- names(event) <- as.character(unique(id))
    pi_u_t <- pi_u_t[names(Time)]

    # subjects who had the event before Thoriz
    ind1 <- Time < Thoriz & event == 1
    # subjects who had the event after Thoriz
    ind2 <- Time > Thoriz
    # subjects who were censored in the interval (Tstart, Thoriz)
    ind3 <- Time < Thoriz & event == 0
    if (any(ind3)) {
        nams <- names(ind3[ind3])
        preds2 <- predict(object, newdata = newdata[id %in% nams, ],
                          process = "event", times = Thoriz, ...)
        weights <- preds2$pred
        f <- factor(preds2$id, levels = unique(preds2$id))
        names(weights) <- f
        weights <- tapply(weights, f, tail, 1)
    }
    loss <- function (x) x * x
    events <- sum(loss(1 - pi_u_t[ind1]), na.rm = TRUE)
    no_events <- sum(loss(pi_u_t[ind2]), na.rm = TRUE)
    censored <- if (any(ind3)) {
        sum(weights * loss(1 - pi_u_t[ind3]) +
                (1 - weights) * loss(pi_u_t[ind3]), na.rm = TRUE)
    } else 0.0
    nr <- length(Time)
    Brier <- (events + no_events + censored) / nr
    out <- list(Brier = Brier, nr = nr, Tstart = Tstart, Thoriz = Thoriz,
                classObject = class(object),
                nameObject = deparse(substitute(object)))
    class(out) <- "tvBrier"
    out
}

print.tvBrier <- function (x, digits = 4, ...) {
    if (!inherits(x, "tvBrier"))
        stop("Use only with 'tvBrier' objects.\n")
    if (x$classObject == "jm")
        cat("\nPrediction Error for the Joint Model", x$nameObject)
    else
        cat("\nPrediction Error for the Cox model", x$nameObject)
    cat("\n\nEstimated Brier score:", round(x$Brier, digits))
    cat("\nAt time:", round(x$Thoriz, digits))
    cat("\nUsing information up to time: ", round(x$Tstart, digits),
        " (", x$nr, " subjects still at risk)", sep = "")
    cat("\n\n")
    invisible(x)
}

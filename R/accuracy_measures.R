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
    if (type_censoring != "right")
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
    if (is.null(newdata[[Time_var]]))
        stop("cannot find the '", Time_var, "' variable in newdata.", sep = "")
    if (is.null(newdata[[event_var]]))
        stop("cannot find the '", event_var, "' variable in newdata.", sep = "")
    newdata <- newdata[newdata[[Time_var]] > Tstart, ]
    newdata <- newdata[newdata[[time_var]] <= Tstart, ]
    newdata[[id_var]] <- newdata[[id_var]][, drop = TRUE]
    test1 <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
    if (!any(test1))
        stop("it seems that there are no events in the interval [Tstart, Thoriz).")
    test2 <- newdata[[Time_var]] > Thoriz & newdata[[event_var]] == 1
    if (!any(test2))
        stop("it seems that there are no events after Thoriz.")
    newdata2 <- newdata
    newdata2[[Time_var]] <- Tstart
    newdata2[[event_var]] <- 0
    preds <- predict(object, newdata = newdata2, process = "event",
                     times = Thoriz, ...)
    pi_u_t <- 1 - preds$pred
    names(pi_u_t) <- preds$id
    pi_u_t <- pi_u_t[preds$times > Tstart]

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
        pi_u_t2 <- preds2$pred
        f <- factor(preds2$id, levels = unique(preds2$id))
        names(pi_u_t2) <- f
        pi_u_t2 <- tapply(pi_u_t2, f, tail, 1)
        nams2 <- names(ind2[ind2])
        ind[ind2] <- ind[ind2] * pi_u_t2[nams2]
    }
    # calculate sensitivity and specificity
    thrs <- seq(0, 1, length = 101)
    Check <- outer(pi_u_t, thrs, "<")
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
    d <- data.frame("cut-off" = x$thrs, "SN" = x$TP, "SP" = 1 - x$FP, "qSN" = x$qSN,
                    "qSP" = x$qSP, "qOverall" = x$qOverall, check.names = FALSE,
                    check.rows = FALSE)
    xx <- rep("", nrow(d))
    xx[ind <- which.max(d$qOverall)] <- "*"
    d[[" "]] <- xx
    nms <- c(as.character(seq(6, 96, 5)), row.names(d)[ind])
    d <- d[row.names(d) %in% nms, ]
    d <- d[!is.na(d$qOverall), ]
    row.names(d) <- 1:nrow(d)
    print(d)
    cat("\n")
    invisible(x)
}

plot.tvROC <- function (x, legend = FALSE,
                        optimal.cutoff = c("", "qOverall", "Youden"), ...) {
    plot(x$FP, x$TP, type = "l", xlab = "1 - Specificity", ylab = "Sensitivity")
    abline(a = 0, b = 1, lty = 3)
    optimal.cutoff <- match.arg(optimal.cutoff)
    if (optimal.cutoff == "qOverall")
        abline(v = x$thrs[which.max(x$qOverall)], lty = 3, lwd = 2, col = 2)
    if (optimal.cutoff == "Youden")
        abline(v = x$thrs[which.max(x$TP - x$FP)], lty = 3, lwd = 2, col = 2)
    if (legend) {
        legend("bottomright", c(paste("At time:", round(x$Thoriz, 1), "\n"),
                                paste("\nUsing information up to time:",
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
    if (type_censoring != "right")
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
    if (is.null(newdata[[Time_var]]))
        stop("cannot find the '", Time_var, "' variable in newdata.", sep = "")
    if (is.null(newdata[[event_var]]))
        stop("cannot find the '", event_var, "' variable in newdata.", sep = "")
    newdata <- newdata[order(newdata[[Time_var]]), ]
    newdata <- newdata[newdata[[Time_var]] > Tstart, ]
    newdata <- newdata[newdata[[time_var]] <= Tstart, ]
    newdata[[id_var]] <- newdata[[id_var]][, drop = TRUE]
    test1 <- newdata[[Time_var]] < Thoriz & newdata[[event_var]] == 1
    if (!any(test1))
        stop("it seems that there are no events in the interval [Tstart, Thoriz).")
    test2 <- newdata[[Time_var]] > Thoriz & newdata[[event_var]] == 1
    if (!any(test2))
        stop("it seems that there are no events after Thoriz.")
    newdata2 <- newdata
    newdata2[[Time_var]] <- Tstart
    newdata2[[event_var]] <- 0
    preds <- predict(object, newdata = newdata2, process = "event",
                     times = Thoriz)
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
                              process = "event", times = Thoriz)
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
                                process = "event", times = Thoriz)
            pi_u_t <- preds4_i$pred
            f <- factor(preds4_i$id, levels = unique(preds4_i$id))
            names(pi_u_t) <- f
            pi_u_t <- tapply(pi_u_t, f, tail, 1)

            preds4_j <- predict(object, newdata = newdata[id %in% unq_nams_j, ],
                                process = "event", times = Thoriz)
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


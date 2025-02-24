tvROC.coxph <-
    function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
              type_weights = c("model-based", "IPCW"), ...) {
    if (!inherits(object, "coxph"))
        stop("Use only with 'coxph' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (!is.null(Thoriz) && Thoriz <= Tstart)
        stop("'Thoriz' must be larger than 'Tstart'.")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    type_censoring <- attr(object$y, "type")
    if (type_censoring != "right")
        stop("'tvROC()' currently only works for right censored data.")
    type_weights <- match.arg(type_weights)
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
                     times = Thoriz, return_mcmc = TRUE, ...)
    qi_u_t <- 1 - preds$mcmc
    rownames(qi_u_t) <- preds$id
    qi_u_t <- qi_u_t[preds$times > Tstart, , drop = FALSE]
    id <- newdata[[id_var]]
    Time <- newdata[[Time_var]]
    event <- newdata[[event_var]]
    f <- factor(id, levels = unique(id))
    Time <- tapply(Time, f, tail, 1L)
    event <- tapply(event, f, tail, 1L)
    names(Time) <- names(event) <- as.character(unique(id))
    thrs <- seq(0, 1, length = 101)
    Check <- lapply(seq_len(ncol(qi_u_t)), function (i) outer(qi_u_t[, i], thrs, "<"))
    Check_mean <- outer(1 - preds$pred[preds$times > Tstart], thrs, "<")
    # outer(qi_u_t, thrs, "<")
    if (type_weights == "model-based") {
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
        ntp <- lapply(Check, function (x) colSums(x * c(ind)))
        tp <- lapply(ntp, function (x) x / sum(ind))
        #nTP <- rowMeans(do.call("cbind", ntp))
        nTP <- colSums(Check_mean * c(ind))
        nFN <- sum(ind) - nTP
        TP <- nTP / sum(ind)
        nfp <- lapply(Check, function (x) colSums(x * c(1 - ind)))
        fp <- lapply(nfp, function (x) x / sum(1 - ind))
        nFP <- rowMeans(do.call("cbind", nfp))
        nFP <- colSums(Check_mean * c(1 - ind))
        nTN <- sum(1 - ind) - nFP
        FP <- nFP / sum(1 - ind)
    } else {
        ind1 <- Time < Thoriz & event == 1
        ind2 <- Time > Thoriz
        nfp <- lapply(Check, function (x) colSums(x * c(ind2)))
        fp <- lapply(nfp, function (x) x / sum(ind2))
        #nFP <- rowMeans(do.call("cbind", nfp))
        nFP <- colSums(Check_mean * c(ind2))
        nTN <- sum(ind2) - nFP
        FP <- nFP / sum(ind2)
        cens_data <- data.frame(Time = Time, cens_ind = 1 - event)
        censoring_dist <- survfit(Surv(Time, cens_ind) ~ 1, data = cens_data)
        weights <- numeric(length(Time))
        ss <- summary(censoring_dist, times = Time[ind1])
        weights[ind1] <- 1 / ss$surv[match(ss$time, Time[ind1])]
        ntp <- lapply(Check, function (x) colSums(x[ind1, , drop = FALSE] / weights[ind1]))
        tp <- lapply(ntp, function (x) x /  sum(1 / weights[ind1]))
        #nTP <- rowMeans(do.call("cbind", ntp))
        nTP <- colSums(Check_mean[ind1, , drop = FALSE] / weights[ind1])
        nFN <- sum(1 / weights[ind1]) - nTP
        TP <- nTP / sum(1 / weights[ind1])
    }
    f1score <- 2 * nTP / (2 * nTP + nFN + nFP)
    F1score <- median(thrs[f1score == max(f1score)])
    youden <- TP - FP
    Youden <- median(thrs[youden == max(youden)])
    out <- list(TP = TP, FP = FP, nTP = nTP, nFN = nFN, nFP = nFP, nTN = nTN,
                tp = do.call("cbind", tp), fp = do.call("cbind", fp),
                thrs = thrs, F1score = F1score, Youden = Youden,
                Tstart = Tstart, Thoriz = Thoriz, nr = length(unique(id)),
                classObject = class(object), type_weights = type_weights,
                nameObject = deparse(substitute(object)))
    class(out) <- "tvROC"
    out
}



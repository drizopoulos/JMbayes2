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
    ind_strata <- attr(object$terms, "specials")$strata 
    LT_MS <- (type_censoring == "counting") & is.null(ind_strata)  
    if (!(type_censoring %in% c("right", "counting")) | !is.null(ind_strata)) 
      stop("'tvROC()' currently only works for (possibly left truncated) right censored data. They are not currently compatible with stratified Cox models.") 
    newdata <- as.data.frame(newdata)
    type_weights <- match.arg(type_weights)
    Tstart <- Tstart + 1e-06
    Thoriz <- Thoriz + 1e-06
    mframe <- model.frame(object$terms, data = newdata)
    Y <- model.response(mframe)
    Time <- if(type_censoring == "right") Y[, 'time'] else Y[, 'stop'] 
    Time0 <- if(type_censoring == "right") rep(0, length(Time)) else Y[, 'start'] 
    newdata <- newdata[Time >= Tstart & Time0 <= Tstart, ]
    mframe <- model.frame(object$terms, data = newdata)
    Y <- model.response(mframe)
    Time <- if(type_censoring == "right") Y[, 'time'] else Y[, 'stop'] 
    Time0 <- if(type_censoring == "right") rep(0, length(Time)) else Y[, 'start'] 
    event <- Y[, 'status']
    sfit0 <- summary(survfit(object, newdata = newdata), times = Tstart)
    sfit <- summary(survfit(object, newdata = newdata), times = Thoriz)
    qi_u_t <- as.matrix(sfit$surv)[1L, ] / as.matrix(sfit0$surv)[1L, ]
    thrs <- seq(0, 1, length = 101)
    Check <- outer(qi_u_t, thrs, "<")
    if (type_weights == "model-based") {
        # subjects who died before Thoriz
        ind1 <- Time < Thoriz & event == 1
        # subjects who were censored in the interval (Tstart, Thoriz)
        ind2 <- Time < Thoriz & event == 0
        ind <- ind1 | ind2
        if (any(ind2)) {
            sfit_cens <- summary(survfit(object, newdata = newdata[ind2, ]),
                             times = Time[ind2])
            sfit2 <- summary(survfit(object, newdata = newdata[ind2, ]),
                             times = Thoriz)
            pi_u_t <- 1.0 - (as.matrix(sfit2$surv)[1L, ] / as.matrix(sfit_cens$surv)[1L, ])
            nams <- names(ind2[ind2])
            ind[nams] <- ind[nams] * pi_u_t[nams]
        }
        # calculate sensitivity and specificity
        nTP <- colSums(Check * c(ind))
        nFN <- sum(ind) - nTP
        TP <- nTP / sum(ind)
        nFP <- colSums(Check * c(1 - ind))
        nTN <- sum(1 - ind) - nFP
        FP <- nFP / sum(1 - ind)
    } else {
        ind1 <- Time < Thoriz & event == 1
        ind2 <- Time > Thoriz
        nFP <- colSums(Check * c(ind2))
        nTN <- sum(ind2) - nFP
        FP <- nFP / sum(ind2)
        cens_data <- data.frame(Time = Time, cens_ind = 1 - event)
        censoring_dist <- survfit(Surv(Time, cens_ind) ~ 1, data = cens_data)
        weights <- numeric(length(Time))
        ss <- summary(censoring_dist, times = Time[ind1])
        weights[ind1] <- 1 / ss$surv[match(ss$time, Time[ind1])]
        nTP <- colSums(Check[ind1, , drop = FALSE] / weights[ind1])
        nFN <- sum(1 / weights[ind1]) - nTP
        TP <- nTP / sum(1 / weights[ind1])
    }
    f1score <- 2 * nTP / (2 * nTP + nFN + nFP)
    F1score <- median(thrs[f1score == max(f1score)])
    youden <- TP - FP
    Youden <- median(thrs[youden == max(youden)])
    out <- list(TP = TP, FP = FP, nTP = nTP, nFN = nFN, nFP = nFP, nTN = nTN,
                thrs = thrs, F1score = F1score, Youden = Youden,
                Tstart = Tstart, Thoriz = Thoriz, nr = nrow(newdata), LT_MS = LT_MS,
                classObject = class(object), type_weights = type_weights,
                nameObject = deparse(substitute(object)))
    class(out) <- "tvROC"
    out
}

tvAUC.coxph <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
                         type_weights = c("model-based", "IPCW"), ...) {
    roc <- tvROC(object, newdata, Tstart, Thoriz, Dt, type_weights, ...)
    TP <- roc$TP
    FP <- roc$FP
    auc <- sum(0.5 * diff(FP) * (TP[-1L] + TP[-length(TP)]), na.rm = TRUE)
    out <- list(auc = auc,
                Tstart = Tstart, Thoriz = roc$Thoriz, nr = roc$nr,
                type_weights = roc$type_weights,
                classObject = class(object),
                LT_MS = roc$LT_MS,
                nameObject = deparse(substitute(object)))
    class(out) <- "tvAUC"
    out
}


tvBrier.coxph <- function (object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
                           integrated = FALSE,
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
    newdata <- as.data.frame(newdata)
    type_censoring <- attr(object$y, "type")
    ind_strata <- attr(object$terms, "specials")$strata 
    LT_MS <- (type_censoring == "counting") & is.null(ind_strata)  
    if (!(type_censoring %in% c("right", "counting")) | !is.null(ind_strata)) 
      stop("'tvBrier()' currently only works for (possibly left truncated) right censored data. They are not currently compatible with stratified Cox models.")
    type_weights <- match.arg(type_weights)
    Tstart <- Tstart + 1e-06
    Thoriz <- Thoriz + 1e-06
    mframe <- model.frame(object$terms, data = newdata)
    Y <- model.response(mframe)
    Time <- if(type_censoring == "right") Y[, 'time'] else Y[, 'stop'] 
    Time0 <- if(type_censoring == "right") rep(0, length(Time)) else Y[, 'start'] 
    newdata <- newdata[Time >= Tstart & Time0 <= Tstart, ]
    mframe <- model.frame(object$terms, data = newdata)
    Y <- model.response(mframe)
    Time <- if(type_censoring == "right") Y[, 'time'] else Y[, 'stop'] 
    Time0 <- if(type_censoring == "right") rep(0, length(Time)) else Y[, 'start'] 
    event <- Y[, 'status']
    brier_fun <- function (pi_u_t, type_weights, weights, ind1, ind2, ind3) {
        loss <- function (x) x * x
        res <- if (type_weights == "model-based") {
            events <- sum(loss(1.0 - pi_u_t[ind1]), na.rm = TRUE)
            no_events <- sum(loss(pi_u_t[ind2]), na.rm = TRUE)
            censored <- if (any(ind3)) {
                sum(weights * loss(1.0 - pi_u_t[ind3]) +
                        (1.0 - weights) * loss(pi_u_t[ind3]), na.rm = TRUE)
            } else 0.0
            (events + no_events + censored) / length(ind1)
        } else {
            mean(loss(as.numeric(ind1) - pi_u_t) * weights)
        }
        res
    }
    ############################################################################
    br <- function (Thoriz) {
        sfit0 <- summary(survfit(object, newdata = newdata), times = Tstart)
        sfit <- summary(survfit(object, newdata = newdata), times = Thoriz)
        pi_u_t <- 1.0 - (as.matrix(sfit$surv)[1L, ] / as.matrix(sfit0$surv)[1L, ])
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
            ss <- summary(censoring_dist, times = Time[ind1])
            weights[ind1] <- 1.0 / ss$surv[match(ss$time, Time[ind1])]
            weights[ind2] <- 1.0 / summary(censoring_dist, times = Thoriz)$surv
        } else {
            weights <- if (any(ind3)) {
                sfit2 <- summary(survfit(object, newdata = newdata[ind3, ]),
                                 times = Thoriz)
                sfit_cens <- summary(survfit(object, newdata = newdata[ind3, ]),
                                 times = Time[ind3])
                1.0 - (as.matrix(sfit2$surv)[1L, ] / as.matrix(sfit_cens$surv)[1L, ] )
            }
        }
        brier_fun(pi_u_t, type_weights, weights, ind1, ind2, ind3)
    }
    Brier <- if (integrated) {
        br1 <- br(0.5 * (Tstart + Thoriz))
        br2 <- br(Thoriz)
        2 * br1 / 3 + br2 / 6
    } else {
        br(Thoriz)
    }
    out <- list(Brier = Brier,
                nr = length(Time), nint = sum(Time < Thoriz & event),
                ncens = sum(Time < Thoriz & event == 0),
                Tstart = Tstart, Thoriz = Thoriz,
                integrated = integrated, type_weights = type_weights, LT_MS = LT_MS,
                classObject = class(object),
                nameObject = deparse(substitute(object)))
    class(out) <- "tvBrier"
    out
}

calibration_plot.coxph <-
    function (object, newdata, Tstart, Thoriz = NULL,
              Dt = NULL, df_ns = NULL, plot = TRUE,
              col = "red", lty = 1, lwd = 1,
              add_CI = TRUE, col_CI = "lightgrey",
              add_density = TRUE, col_dens = "grey",
              xlab = "Predicted Probabilities",
              ylab = "Observed Probabilities", main = "", ...) {
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
        stop("'newdata' must be a data.frame with more than one rows.\n")
    if (is.null(Thoriz) && is.null(Dt))
        stop("either 'Thoriz' or 'Dt' must be non null.\n")
    if (!is.null(Thoriz) && Thoriz <= Tstart)
        stop("'Thoriz' must be larger than 'Tstart'.")
    if (is.null(Thoriz))
        Thoriz <- Tstart + Dt
    newdata <- as.data.frame(newdata)
    type_censoring <- object$model_info$type_censoring
    Tstart <- Tstart + 1e-06
    Thoriz <- Thoriz + 1e-06
    mframe <- model.frame(object$terms, data = newdata)
    Y <- model.response(mframe)
    Time <- Y[, 'time']
    newdata <- newdata[Time >= Tstart, ]
    mframe <- model.frame(object$terms, data = newdata)
    Y <- model.response(mframe)
    Time <- Y[, 'time']
    event <- Y[, 'status']
    sfit0 <- summary(survfit(object, newdata = newdata), times = Tstart)
    sfit <- summary(survfit(object, newdata = newdata), times = Thoriz)
    pi_u_t <- 1.0 - (as.matrix(sfit$surv)[1L, ] / as.matrix(sfit0$surv)[1L, ])
    cloglog <- function (x) log(-log(1.0 - x))
    cal_DF <- data.frame(Time = pmin(Thoriz, Time),
                         event = event * (Time <= Thoriz), preds = pi_u_t)
    if (is.null(df_ns)) {
        df_ns <- if (sum(cal_DF$event) > 25L) 3L else 2L
    } else {
        if (!is.numeric(df_ns) || length(df_ns) != 1L)
            stop("'df_ns' must be a numeric scalar.")
    }
    form <- paste0("Surv(Time, event) ~ nsk(cloglog(preds), df = ", df_ns, ")")
    cal_Cox <- coxph(as.formula(form), data = cal_DF)
    qs <- quantile(pi_u_t, probs = c(0.01, 0.99))
    probs_grid <- data.frame(preds = seq(qs[1L], qs[2L], length.out = 101L))
    ss <- summary(survfit(cal_Cox, newdata = probs_grid), times = Thoriz)
    obs <- 1 - c(ss$surv)
    low <- 1 - c(ss$low)
    upp <- 1 - c(ss$upp)
    obs_pi_u_t <- 1 - c(summary(survfit(cal_Cox, newdata = cal_DF),
                                times = Thoriz)$surv)
    preds <- probs_grid$preds
    if (plot) {
        plot(preds, obs, type = "n", xlab = xlab, ylab = ylab, main = main,
             xlim = c(0, 1), ylim = c(0, 1))
        if (add_CI) {
            polygon(c(preds, rev(preds)), c(low, rev(upp)), border = NA,
                    col = col_CI)
        }
        lines(preds, obs, lty = lty, col = col, lwd = lwd)
        abline(0, 1, lty = 2)
        if (add_density) {
            par(new = TRUE)
            plot(density(pi_u_t), axes = FALSE, xlab = "", ylab = "", main = "",
                 col = col_dens)
        }
        invisible()
    } else {
        list("observed" = obs, "predicted" = preds, "pi_u_t" = pi_u_t,
             "obs_pi_u_t" = obs_pi_u_t, "low" = low, "upp" = upp)
    }
}



#############################################################################
# Clinical utility metrics and decision curves
tvMetric.coxph <- function (object, newdata, Tstart, 
                            metric = c("IWEC", "NB", "EC"),
                            shape1 = NULL, shape2 = NULL, mod = NULL, interv = NULL, mass = 0.95,
                            decision_threshold = NULL, Thoriz = NULL, Dt = NULL,
                            type_weights = c("model-based", "IPCW"), 
                            opt_interv = c(1, 500), name_model = NULL, ...) {
  if (!inherits(object, "coxph"))
    stop("Use only with 'coxph' objects.\n")
  if (!is.data.frame(newdata) || nrow(newdata) == 0)
    stop("'newdata' must be a data.frame with more than one rows.\n")
  if (is.null(Thoriz) && is.null(Dt))
    stop("either 'Thoriz' or 'Dt' must be non null.\n")
  if (!is.null(Thoriz) && Thoriz <= Tstart)
    stop("'Thoriz' must be larger than 'Tstart'.")
  if(metric == "IWEC") {
    if ( (is.null(shape1) | is.null(shape2)) & (is.null(mod) | is.null(interv) | is.null(mass) ) ) {
      stop("either 'shape1' and 'shape2', or 'mod', 'interv' and 'mass' must be non null.\n")
    } 
    if ( !is.null(shape1) & !is.null(shape2) ) {
      if ( shape1 < 0 || shape2 < 0 )
        stop("shape1 and shape2 must be greater or equal to 0.\n")
    }
    if ( !is.null(mod) & !is.null(interv) & !is.null(mass) ) {
      if ( mod < 0 | mod > 1 )
        stop("'mod' must be between 0 and 1.\n")
      if ( length(interv) != 2 ) {
        stop("'interv' must be a vector with two elements.\n")
      } else if( interv[2] < interv[1] ) {
        stop("The second element of 'interv' must be greater than the first.\n")
      } else if( (interv[1] < 0 | interv[1] > 1) |  (interv[2] < 0 | interv[2] > 1)) {
        stop("The two element in 'interv' must be between 0 and 1.\n")
      }
      if ( mass < 0 | mass > 1 )
        stop("'mass' must be between 0 and 1.")
    } 
    w_function_builder <- function(mod, ext1, ext2, p, interval = opt_interv) {
      f <- function(x) {
        y <- (x - 1) / mod - x + 2
        pbeta(ext2, x, y) - pbeta(ext1, x, y) - p
      }
      root <- try(uniroot(f, interval = interval)$root, TRUE)
      if (inherits(root, "try-error"))
        stop("Failed to compute shape parameters because no valid solution could be found. Please check 'mod', 'interv', and 'mass'.\n")
      y <- (root - 1) / mod - root + 2
      c(root, y)
    }
    if ( is.null(shape1) & is.null(shape2) ) {
      shapes <- w_function_builder(mod, interv[1], interv[2], mass)
      shape1 <- shapes[1]
      shape2 <- shapes[2]
    }
  } else {
    if(is.null(decision_threshold))
      stop("When 'NB' or 'EC' are used, 'decision_threshold' must non-null.")
    if (decision_threshold <= 0 || decision_threshold >= 1 )
      stop("decision_threshold must be between 0 and 1.")
  }
  if (is.null(Thoriz))
    Thoriz <- Tstart + Dt
  newdata <- as.data.frame(newdata)
  type_censoring <- attr(object$y, "type")
  ind_strata <- attr(object$terms, "specials")$strata 
  LT_MS <- (type_censoring == "counting") & is.null(ind_strata)  
  if (!(type_censoring %in% c("right", "counting")) | !is.null(ind_strata)) 
    stop("Predictive performance metrics currently only work for (possibly left truncated) right censored data. They are not currently compatible with stratified Cox models.") 
  type_weights <- match.arg(type_weights)
  Tstart <- Tstart + 1e-06
  Thoriz <- Thoriz + 1e-06
  mframe <- model.frame(object$terms, data = newdata)
  Y <- model.response(mframe)
  Time <- if(type_censoring == "right") Y[, 'time'] else Y[, 'stop'] 
  Time0 <- if(type_censoring == "right") rep(0, length(Time)) else Y[, 'start'] 
  newdata <- newdata[Time >= Tstart & Time0 <= Tstart, ]
  mframe <- model.frame(object$terms, data = newdata)
  Y <- model.response(mframe)
  Time <- if(type_censoring == "right") Y[, 'time'] else Y[, 'stop'] 
  Time0 <- if(type_censoring == "right") rep(0, length(Time)) else Y[, 'start'] 
  event <- Y[, 'status']
  sfit0 <- summary(survfit(object, newdata = newdata), times = Tstart)
  sfit <- summary(survfit(object, newdata = newdata), times = Thoriz)
  pi_u_t <- 1.0 - (as.matrix(sfit$surv)[1L, ] / as.matrix(sfit0$surv)[1L, ])
  w_MB <- NULL
  if(type_weights == "model-based") {
    ind3 <- Time <= Thoriz & event == 0
    w_MB <- rep(NA, length(Time))
    if (any(ind3)) {
      sfit2 <- summary(survfit(object, newdata = newdata[ind3, ]),
                       times = Thoriz)
      sfit_cens <- summary(survfit(object, newdata = newdata[ind3, ]),
                           times = Time[ind3])
      w_MB[ind3] <- 1.0 - (as.matrix(sfit2$surv)[1L, ] / as.matrix(sfit_cens$surv)[1L, ])
    }
  } 
  out <- tvMetric.default(pi_u_t, Time, event, LT_Time = Time0, Tstart = Tstart, 
                          metric = metric, shape1 = shape1, shape2 = shape2, 
                          mod = mod, interv = interv, mass = mass,
                          decision_threshold = decision_threshold, Thoriz = Thoriz, 
                          type_weights = type_weights, preds_MB = w_MB, 
                          opt_interv = opt_interv)
  metric_info <- if(metric == "IWEC") {
    c(shape1, shape2)
  } else {
    decision_threshold
  }
  if(is.null(name_model)) {
    nameObject <- deparse(substitute(object))
  } else {
    if(!is.character(name_model))
      stop("'name_model' must be a character vector")
    nameObject <- name_model
  }
  out <- list(metric = out$metric,
              metric_type = metric, 
              nr = out$nr, nint = out$nint,
              ncens = out$ncens,
              Tstart = Tstart, Thoriz = Thoriz,
              type_weights = type_weights,
              metric_info = metric_info,
              LT_MS = LT_MS, 
              classObject = class(object),
              nameObject = nameObject)
  class(out) <- "tvMetric"
  out
}

## Metric wrappers:
tvIWEC.coxph <- function(object, newdata, Tstart, 
                         shape1 = NULL, shape2 = NULL, mod = NULL, interv = NULL, mass = 0.95,
                         Thoriz = NULL, Dt = NULL,
                         type_weights = c("model-based", "IPCW"), 
                         opt_interv = c(1, 500), name_model = NULL, ...) {
  name_model <- if(is.null(name_model)) {
    deparse(substitute(object))
  } else {
    name_model
  } 
  tvMetric.coxph(object, newdata = newdata, Tstart = Tstart, 
                 metric = "IWEC", shape1 = shape1, shape2 = shape2, 
                 mod = mod, interv = interv, mass = mass,
                 decision_threshold = NULL, Thoriz = Thoriz, Dt = Dt,
                 type_weights = type_weights, opt_interv = opt_interv, name_model = name_model, ...)
}

tvNB.coxph <- function(object, newdata, Tstart, 
                       decision_threshold = NULL, Thoriz = NULL, Dt = NULL,
                       type_weights = c("model-based", "IPCW"), name_model = NULL, ...) {
  name_model <- if(is.null(name_model)) {
    deparse(substitute(object))
  } else {
    name_model
  } 
  tvMetric.coxph(object = object, newdata = newdata, Tstart = Tstart, 
                 metric = "NB", shape1 = NULL, shape2 = NULL, 
                 mod = NULL, interv = NULL, mass = NULL,
                 decision_threshold = decision_threshold, Thoriz = Thoriz, Dt = Dt,
                 type_weights = type_weights, opt_interv = NULL, name_model = name_model, ...)
}

tvEC.coxph <- function(object, newdata, Tstart, 
                       decision_threshold = NULL, Thoriz = NULL, Dt = NULL,
                       type_weights = c("model-based", "IPCW"), name_model = NULL, ...) {
  name_model <- if(is.null(name_model)) {
    deparse(substitute(object))
  } else {
    name_model
  } 
  tvMetric.coxph(object = object, newdata = newdata, Tstart = Tstart, 
                 metric = "EC", shape1 = NULL, shape2 = NULL, 
                 mod = NULL, interv = NULL, mass = NULL,
                 decision_threshold = decision_threshold, Thoriz = Thoriz, Dt = Dt,
                 type_weights = type_weights, opt_interv = NULL, name_model = name_model, ...)
}


decision_curves.coxph <-
  function (object, newdata, Tstart, 
            decision_threshold_range = seq(0 + 1e-06, 1 - 1e-06, length.out = 51),
            type = c("EC", "weightedEC", "NB"),
            shape1 = NULL, shape2 = NULL, mod = NULL, interv = NULL, mass = 0.95, 
            Thoriz = NULL, Dt = NULL, type_weights = c("model-based", "IPCW"),
            opt_interv = c(1, 500), ...) {
    if (!inherits(object, "coxph"))
      stop("Use only with 'coxph' objects.\n")
    if (!is.data.frame(newdata) || nrow(newdata) == 0)
      stop("'newdata' must be a data.frame with more than one rows.\n")
    newdata_original <- newdata
    if (is.null(Thoriz) && is.null(Dt)) {
      stop("either 'Thoriz' or 'Dt' must be non null.\n")
    }
    if (!is.null(Thoriz) && Thoriz <= Tstart) {
      stop("'Thoriz' must be larger than 'Tstart'.")
    }
    if (length(decision_threshold_range) < 2 || any(decision_threshold_range <= 0) || any(decision_threshold_range >= 1))
      stop("Multiple thresholds are required and must be between 0 and 1.")
    if (is.null(Thoriz)) {
      Thoriz <- Tstart + Dt
    }
    type <- match.arg(type)
    if ((type == "weightedEC") && ((is.null(shape1) || is.null(shape2)) && (is.null(mod) || is.null(interv) || is.null(mass) )) ) {
      stop("If type is 'weightedEC', either 'shape1' and 'shape2', or 'mod', 'interv' and 'mass' must be non null.\n")
    } 
    if ( (type == "weightedEC") && (!is.null(shape1) && !is.null(shape2)) ) {
      if ( shape1 < 0 || shape2 < 0 )
        stop("shape1 and shape2 must be greater or equal to 0.\n")
    }
    if ((type == "weightedEC") & (!is.null(mod) & !is.null(interv) & !is.null(mass)) ) {
      if ( mod < 0 | mod > 1 )
        stop("'mod' must be between 0 and 1.\n")
      if ( length(interv) != 2 ) {
        stop("'interv' must be a vector with two elements.\n")
      } else if( interv[2] < interv[1] ) {
        stop("The second element of 'interv' must be greater than the first.\n")
      } else if( (interv[1] < 0 | interv[1] > 1) |  (interv[2] < 0 | interv[2] > 1)) {
        stop("The two element in 'interv' must be between 0 and 1.\n")
      }
      if ( mass < 0 | mass > 1 )
        stop("'mass' must be between 0 and 1.")
    } 
    w_function_builder <- function(mod, ext1, ext2, p, interval = opt_interv) {
      f <- function(x) {
        y <- (x - 1)/mod - x + 2
        pbeta(ext2, x, y) - pbeta(ext1, x, y) - p
      }
      root <- try(uniroot(f, interval = interval)$root, TRUE)
      if (inherits(root, "try-error"))
        stop("Failed to compute shape parameters because no valid solution could be found. Please check 'mod', 'interv', and 'mass'.\n")
      y <- (root - 1)/mod - root + 2
      c(root, y)
    }
    if ( (type == "weightedEC") & (is.null(shape1) & is.null(shape2) )) {
      shapes <- w_function_builder(mod, interv[1], interv[2], mass)
      shape1 <- shapes[1]
      shape2 <- shapes[2]
    }
    newdata <- as.data.frame(newdata)
    type_censoring <- attr(object$y, "type")
    ind_strata <- attr(object$terms, "specials")$strata #!
    ## with the current implementation left-truncated and right-censored data cannot
    ## be used with stratified Cox model
    LT_MS <- (type_censoring == "counting") & is.null(ind_strata) #! 
    if (!(type_censoring %in% c("right", "counting")) | !is.null(ind_strata)) #!
      stop("'decision_curves()' currently only works for (possibly left truncated) right censored data. It is not currently compatible with stratified Cox models.") #!
    type_weights <- match.arg(type_weights)
    Tstart <- Tstart + 1e-06
    Thoriz <- Thoriz + 1e-06
    mframe <- model.frame(object$terms, data = newdata)
    Y <- model.response(mframe)
    Time <- if(type_censoring == "right") Y[, 'time'] else Y[, 'stop'] #!
    Time0 <- if(type_censoring == "right") rep(0, length(Time)) else Y[, 'start'] #!
    newdata <- newdata[Time >= Tstart & Time0 <= Tstart, ]
    mframe <- model.frame(object$terms, data = newdata)
    Y <- model.response(mframe)
    Time <- if(type_censoring == "right") Y[, 'time'] else Y[, 'stop'] #!
    Time0 <- if(type_censoring == "right") rep(0, length(Time)) else Y[, 'start'] #!
    event <- Y[, 'status']
    if(type == "NB"){
      loss <- function(pi, di, c) as.numeric(pi > c) * (di - ((1 - di) * (c)) / (1-c))
      loss_all <- function(di, c) di - (1 - di) * (c / (1 - c))
      loss_none <- function(di, c) 0
    } else {
      loss <- function(pi, di, c) as.numeric(pi > c) * c * (1 - di) + (1 - as.numeric(pi > c)) * (1 - c) * di
      loss_all <- function(di, c) c * (1 - di)
      loss_none <- function(di, c) (1 - c) * di
    }
    lo <- function (Thoriz) {
      sfit0 <- summary(survfit(object, newdata = newdata), times = Tstart)
      sfit <- summary(survfit(object, newdata = newdata), times = Thoriz)
      pi_u_t <- 1.0 - (as.matrix(sfit$surv)[1L, ] / as.matrix(sfit0$surv)[1L, ])
      # number of subjects ar risk
      nr <- length(Time)
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
        cens_data <- data.frame(Time0 = Time0, Time = Time, cens_ind = 1 - event)
        censoring_dist <- survfit(Surv(Time0, Time, cens_ind) ~ 1, data = cens_data)
        weights <- numeric(length(Time))
        ss <- summary(censoring_dist, times = Time[ind1])
        weights[ind1] <- 1.0 / ss$surv[match(ss$time, Time[ind1])]
        weights[ind2] <- 1.0 / summary(censoring_dist, times = Thoriz)$surv
        res <- list(3)
        res[[1]] <- vapply(decision_threshold_range, function(cc) mean(loss(pi_u_t, di = as.numeric(ind1), c = cc) * weights),
                           numeric(1))
        res[[2]] <- vapply(decision_threshold_range, function(cc) mean(loss_all(di = as.numeric(ind1), c = cc) * weights),
                           numeric(1))
        res[[3]] <- vapply(decision_threshold_range, function(cc) mean(loss_none(di = as.numeric(ind1), c = cc) * weights),
                           numeric(1))
      } else {
        preds_MB <- rep(NA, length(Time))
        if (any(ind3)) {
          sfit2 <- summary(survfit(object, newdata = newdata[ind3, ]),
                           times = Thoriz)
          sfit_cens <- summary(survfit(object, newdata = newdata[ind3, ]),
                               times = Time[ind3])
          weights <- 1.0 - (as.matrix(sfit2$surv)[1L, ] / as.matrix(sfit_cens$surv)[1L, ] )
          preds_MB[ind3] <- weights
        }
        res <- list(3)
        events_model <- vapply(decision_threshold_range, function(cc) sum(loss(pi_u_t[ind1], di = as.numeric(ind1)[ind1], 
                                                                               c = cc), na.rm = TRUE),
                               numeric(1))
        events_all <- vapply(decision_threshold_range, function(cc) sum(loss_all(di = as.numeric(ind1)[ind1], 
                                                                                 c = cc), na.rm = TRUE),
                             numeric(1))
        events_none <- vapply(decision_threshold_range, function(cc) sum(loss_none(di = as.numeric(ind1)[ind1], 
                                                                                   c = cc), na.rm = TRUE),
                              numeric(1))
        no_events_model <- vapply(decision_threshold_range, function(cc) sum(loss(pi_u_t[ind2], di = as.numeric(ind1)[ind2], 
                                                                                  c = cc), na.rm = TRUE),
                                  numeric(1))
        no_events_all <- vapply(decision_threshold_range, function(cc) sum(loss_all(di = as.numeric(ind1)[ind2], 
                                                                                    c = cc), na.rm = TRUE),
                                numeric(1))
        no_events_none <- vapply(decision_threshold_range, function(cc) sum(loss_none(di = as.numeric(ind1)[ind2], 
                                                                                      c = cc), na.rm = TRUE),
                                 numeric(1))
        censored_model <- if (any(ind3)) {
          vapply(decision_threshold_range, function(cc) sum(weights * loss(pi_u_t[ind3], di = 1.0, c = cc) +
                                                              (1 - weights) * loss(pi_u_t[ind3], di = 0.0,
                                                                                   c = cc), na.rm = TRUE),
                 numeric(1))
        } else rep(0.0, length(decision_threshold_range))
        censored_all <- if (any(ind3)) {
          vapply(decision_threshold_range, function(cc) sum(weights * loss_all(di = 1.0, c = cc) +
                                                              (1 - weights) * loss_all(di = 0.0,
                                                                                       c = cc), na.rm = TRUE),
                 numeric(1))
        } else rep(0.0, length(decision_threshold_range))
        censored_none <- if (any(ind3)) {
          vapply(decision_threshold_range, function(cc) sum(weights * loss_none(di = 1.0, c = cc) +
                                                              (1 - weights) * loss_none(di = 0.0,
                                                                                        c = cc), na.rm = TRUE),
                 numeric(1))
        } else rep(0.0, length(decision_threshold_range))
        res[[1]] <- (events_model + no_events_model + censored_model) / length(ind1)
        res[[2]] <- (events_all + no_events_all + censored_all) / length(ind1)
        res[[3]] <- (events_none + no_events_none + censored_none) / length(ind1)
      }
      ## compute IWECs to be printed in the plot
      iwecs <- numeric(3)
      if(type == "weightedEC") {
        iwecs[1] <- tvMetric(pi_u_t, Time = Time, event = event, LT_Time = Time0, Tstart = Tstart,
                             metric = "IWEC", shape1 = shape1, shape2 = shape2, Thoriz = Thoriz,
                             type_weights = type_weights, preds_MB = preds_MB, opt_interv)$metric
        iwecs[2] <- tvMetric(rep(1, length(Time)), Time = Time, event = event, LT_Time = Time0, Tstart = Tstart,
                             metric = "IWEC", shape1 = shape1, shape2 = shape2, Thoriz = Thoriz,
                             type_weights = type_weights, preds_MB = preds_MB, opt_interv)$metric
        iwecs[3] <- tvMetric(rep(0, length(Time)), Time = Time, event = event, LT_Time = Time0, Tstart = Tstart,
                             metric = "IWEC", shape1 = shape1, shape2 = shape2, Thoriz = Thoriz,
                             type_weights = type_weights, preds_MB = preds_MB, opt_interv)$metric
      }
      list(metric_c = res, nr = nr, ind1 = ind1, ind2 = ind2, ind3 = ind3, iwecs = iwecs)
    }
    out <- lo(Thoriz)
    out <- list(metric_c = out$metric_c,
                newdata = newdata_original,
                nr = out$nr, nint = sum(out$ind1),
                ncens = sum(out$ind3),
                Tstart = Tstart, Thoriz = Thoriz,
                type = type, shape1 = shape1, shape2 = shape2,
                decision_threshold_range = decision_threshold_range,
                type_weights = type_weights,
                iwecs = out$iwecs,
                classObject = class(object),
                nameObject = deparse(substitute(object)))
    class(out) <- "decision_curves"
    out
  }
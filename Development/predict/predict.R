pred_Long <- predLong2
pred_Event <- predSurv
subject <- 1
outcomes <- 1
CI <- TRUE
xlab <- "time"
ylab_long <- NULL
lwd_long <- 2
lwd_event <- 2
col_line_long <- "blue"
col_line_event <- "red"
fill_CI_long <- "#0000FF80"
fill_CI_event <- "#FF000080"
cex_xlab <- 1
cex_ylab_long <- 1

####

id_var <- "id"
time_var <- "year"
resp_vars <- jointFit$model_info$var_names$respVars_form
ranges <- lapply(jointFit$model_data$y, range, na.rm = TRUE)
last_times <- rep(3, 2)

####

test1 <- is.data.frame(pred_Long)
test2 <- is.list(pred_Long) && length(pred_Long) == 2L && is.data.frame(pred_Long[[1]])
if (!test1 && !test2) {
    stop("'pred_Long' must be the output of ",
         "predict.jm(..., return_newdata = TRUE)")
}
if (test2) {
    pred_Long <- rbind(pred_Long[[1L]], pred_Long[[2L]])
}
unq_id <- unique(pred_Long[[id_var]])
if (length(subject) > 1L) {
    stop("'subject' must be of length 1.")
}
if (!subject %in% unq_id && subject > length(unq_id)) {
    stop("not valid input for 'subject'.")
}
subj <- if (subject %in% unq_id) subject else unq_id[subject]
pred_Long <- pred_Long[pred_Long[[id_var]] == subj, ]
subj_ind <- match(subj, unq_id)
pred_Long <- pred_Long[pred_Long[[time_var]] <= last_times[subj_ind], ]
pred_Event <- pred_Event[pred_Event[[id_var]] == subj, ]
pos_outcomes <- grep("pred_", names(pred_Long), fixed = TRUE)
n_outcomes <- length(pos_outcomes)
if (any(outcomes > n_outcomes)) {
    stop("not valid entries in 'outcome'.")
}
if (is.null(ylab_long)) {
    ylab_long <- resp_vars
}

xlim <- NULL
xlim <- if (!is.null(pred_Long)) range(xlim, pred_Long[[time_var]])
xlim <- if (!is.null(pred_Event)) range(xlim, pred_Event[[time_var]])

plot_long_i <- function (outcome, add_xlab = FALSE, box = TRUE) {
    ind <- pos_outcomes[outcome]
    preds <- pred_Long[[ind]]
    low <- pred_Long[[ind + 1]]
    upp <- pred_Long[[ind + 2]]
    times <- pred_Long[[time_var]]
    ry <- if (CI) range(preds, low, upp) else range(preds)
    rx <- range(times)
    plot(rx, ry, type = "n", xaxt = "n", bty = if (box) "o" else "n",
         xlab = if (add_xlab) xlab  else "", xlim = xlim,
         ylim = ranges[[outcome]], ylab = ylab_long[outcome],
         cex.lab = cex_ylab_long)
    if (CI) {
        polygon(c(times, rev(times)), c(low, rev(upp)), border = NA,
                col = fill_CI_long)
    }
    lines(pred_Long[[time_var]], pred_Long[[ind]],
          lwd = lwd_long, col = col_line_long)
    if (!is.null(pred_Event)) abline(v = last_times[subj_ind] + 0.01, lty = 3)
}
plot_event <- function (box = FALSE) {
  ind <- grep("pred_", names(pred_Event), fixed = TRUE)
  preds <- pred_Event[[ind]]
  low <- pred_Event[[ind + 1]]
  upp <- pred_Event[[ind + 2]]
  times <- pred_Event[[time_var]]
  rx <- range(times)
  plot(rx, c(0, 1), type = "n", xlab = "", ylab = "", xlim = xlim, axes = FALSE)
  if (box) box()
  axis(4)
  if (CI) {
    polygon(c(times, rev(times)), c(low, rev(upp)), border = NA,
            col = fill_CI_event)
  }
  lines(pred_Event[[time_var]], pred_Event[[ind]],
        lwd = lwd_long, col = col_line_event)
}


# n_outcomes == 1
op <- par(mar = c(4,4,4,4), mgp = c(2, 0.4, 0), tcl = -0.3)
plot_long_i(1, TRUE)
axis(1)
par(new = TRUE)
plot_event()
mtext("CIF", 4, 2)
par(op)


# n_outcomes == 2
op <- par(mfrow = c(2, 1), oma = c(4,4,4,4), mar = c(0, 0, 0, 0),
          mgp = c(2, 0.4, 0), tcl = -0.3)
plot_long_i(1, box = FALSE)
axis(1, c(-5, last_times[subj_ind]), labels = c("", ""), tcl = 0)
plot_long_i(2, box = FALSE)
axis(1)
mtext(xlab, side = 1, line = 1.5, outer = TRUE, cex = cex_xlab)
par(op)
op <- par(new = TRUE, oma = c(4,4,4,4), mar = c(0, 0, 0, 0),
          mgp = c(2, 0.4, 0), tcl = -0.3, cex = 0.9)
plot_event(box = TRUE)
mtext("CIF", 4, 2)
par(op)

# n_outcomes == 3
op <- par(mfrow = c(3, 1), oma = c(4,4,4,4), mar = c(0, 0, 0, 0),
          mgp = c(2, 0.4, 0), tcl = -0.3)
pos <- par("usr")[1] + c(0.25, 0.5, 0.75) * diff(par("usr")[1:2])
plot_long_i(1, box = FALSE)
axis(1, c(-5, last_times[subj_ind]), labels = c("", ""), tcl = 0)
mtext("serBilir", 2, 1.7, at = pos[1])
plot_long_i(2, box = FALSE)
axis(1, c(-5, last_times[subj_ind]), labels = c("", ""), tcl = 0)
mtext("Prothro", 2, 1.7, at = pos[2])
plot_long_i(3, box = FALSE)
mtext("Ascites", 2, 1.7, at = pos[3])
axis(1)
mtext(xlab, side = 1, line = 1.5, outer = TRUE, cex = cex_xlab)
box("inner")
par(op)

op <- par(new = TRUE, oma = c(1.9, 2.61, 2, 2.61), mar = c(0, 0, 0, 0),
          mgp = c(2, 0.4, 0), tcl = -0.3, cex = 0.66)
plot_event()
mtext("CIF", 4, 1.5)
par(op)














################################################################################
################################################################################
################################################################################

if (FALSE) {
    library("JMbayes2")
    pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
    CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
    fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
    fm2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)
    fm3 <- mixed_model(ascites ~ year + sex, data = pbc2,
                       random = ~ year | id, family = binomial())

    jointFit1 <- jm(CoxFit, list(fm1), time_var = "year")

    source("./R/help_functions.R")
    source("./R/predict_funs.R")
    Rcpp::sourceCpp('src/mcmc_fit.cpp')
}


object <- jointFit1
ND <- pbc2[pbc2$id %in% c(2, 3, 15), ]
ND$id <- factor(ND$id)
ND2 <- ND[ND$year > 1, ]
ND <- ND[ND$year < 1, ]
ND$status2 <- 0
ND$years <- with(ND, ave(year, id, FUN = function (x) max(x, na.rm = T)))
ND$prothrombin[c(3, 5, 8)] <- NA
newdata = ND
newdata2 = ND2
times = NULL
process = "event"
type_pred = "response"
type = "subject_specific"
level = 0.95; return_newdata = FALSE
n_samples = 500L; n_mcmc = 55L; cores = NULL; seed = 123L

#############################################################
#############################################################

if (is.null(cores)) {
    n <- length(unique(newdata[[object$model_info$var_names$idVar]]))
    cores <- if (n > 20) 4L else 1L
}
components_newdata <- get_components_newdata(object, newdata, n_samples,
                                             n_mcmc, cores, seed)

components_newdata$mcmc$b[, , 21]

object$model_info$CR_MS

predict_Event <- function (object, components_newdata, newdata, level) {
    control <- object$control
    terms_FE <- object$model_info$terms$terms_FE
    terms_FE_noResp <- object$model_info$terms$terms_FE_noResp
    terms_RE <- object$model_info$terms$terms_RE
    idVar <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    terms_Surv <- object$model_info$terms$terms_Surv
    terms_Surv_noResp <- object$model_info$terms$terms_Surv_noResp
    type_censoring <- object$model_info$type_censoring
    dataL <- newdata
    Xbar <- object$model_data$Xbar
    data_pred <- newdata
    idT <- data_pred[[idVar]]
    data_pred <- data_pred[tapply(row.names(data_pred),
                                  factor(idT, unique(idT)), tail, 1L), ]
    mf_data_pred <- model.frame.default(terms_Surv, data = data_pred)
    Surv_Response <- model.response(mf_data_pred)
    ind_strata <- attr(terms_Surv, "specials")$strata
    strata <- if (is.null(ind_strata)) {
        rep(1, nrow(mf_data_pred))
    } else {
        unclass(mf_data_pred[[ind_strata]])
    }
    # The definition of last_times needs to be checked for counting and interval
    last_times <- switch(type_censoring, "right" = unname(Surv_Response[, "time"]),
                         "counting" = unname(Surv_Response[, "stop"]),
                         "interval" = unname(Surv_Response[, "time1"]))
    t_max <- quantile(object$model_data$Time_right, probs = 0.9)
    # times <- seq(0.5, 20, len = 50)
    if (is.null(times) || !is.numeric(times)) {
        times <- lapply(last_times, seq, to = t_max, length.out = 21L)
    } else {
        t_max <- max(object$model_data$Time_right)
        f <- function (lt, tt, tm) c(lt, sort(tt[tt > lt & tt <= tm]))
        times <- lapply(last_times, f, tt = times, tm = t_max)
    }

    n_times <- sapply(times, length)
    data_pred <- data_pred[rep(seq_along(times), n_times), ]
    data_pred[[time_var]] <- unlist(times, use.names = FALSE)
    idT <- data_pred[[idVar]]
    idT <- factor(idT, levels = unique(idT))
    strata <- rep(strata, n_times)
    upp_limit <- data_pred[[time_var]]
    Time_start <- last_times[unclass(idT)]
    g <- function (t0, t) c(t0, head(t, -1))
    low_limit <- unlist(mapply2(g, last_times, times), use.names = FALSE)
    GK <- gaussKronrod(k = 7L)
    sk <- GK$sk
    P <- c(upp_limit - low_limit) / 2
    st <- outer(P, sk) + (c(upp_limit + low_limit) / 2)
    log_Pwk <- unname(rep(log(P), each = length(sk)) +
                          rep_len(log(GK$wk), length.out = length(st)))

    # knots
    knots <- control$knots

    # indices
    ni_event <- tapply(idT, idT, length)
    ni_event <- cbind(c(0, head(cumsum(ni_event), -1)), cumsum(ni_event))
    id_H <- rep(paste0(idT, "_", unlist(tapply(idT, idT, seq_along))),
                each = 7L)
    id_H <- match(id_H, unique(id_H))
    # id_H_ repeats each unique idT the number of quadrature points
    id_H_ <- rep(idT, each = 7L)
    id_H_ <- match(id_H_, unique(id_H_))
    id_h <- unclass(idT)


    # Functional forms
    functional_forms <- object$model_info$functional_forms
    FunForms_per_outcome <- object$model_info$FunForms_per_outcome
    collapsed_functional_forms <- object$model_info$collapsed_functional_forms
    FunForms_cpp <- object$model_info$FunForms_cpp
    FunForms_ind <- object$model_info$FunForms_ind
    Funs_FunForms <- object$model_info$Funs_FunForms
    eps <- object$model_info$eps
    direction <- object$model_info$direction


    strata_H <- rep(strata, each = 7L)
    W0_H <- create_W0(c(t(st)), knots, control$Bsplines_degree + 1, strata_H)
    dataS_H <- SurvData_HazardModel(st, data_pred, Time_start,
                                    paste0(idT, "_", strata), time_var)
    mf <- model.frame.default(terms_Surv_noResp, data = dataS_H)
    W_H <- construct_Wmat(terms_Surv_noResp, mf)
    any_gammas <- as.logical(ncol(W_H))
    if (!any_gammas) {
        W_H <- matrix(0.0, nrow = nrow(W_H), ncol = 1L)
    }
    attr <- lapply(functional_forms, extract_attributes, data = dataS_H)
    eps <- lapply(attr, "[[", 1L)
    direction <- lapply(attr, "[[", 2L)
    X_H <- design_matrices_functional_forms(st, terms_FE_noResp,
                                            dataL, time_var, idVar, idT,
                                            collapsed_functional_forms, Xbar,
                                            eps, direction)
    Z_H <- design_matrices_functional_forms(st, terms_RE,
                                            dataL, time_var, idVar, idT,
                                            collapsed_functional_forms, NULL,
                                            eps, direction)
    U_H <- lapply(functional_forms, construct_Umat, dataS = dataS_H)

    X_H[] <- lapply(X_H, docall_cbind)
    Z_H[] <- lapply(Z_H, docall_cbind)

    Data <- list(
        log_Pwk = log_Pwk, id_H = id_H, id_h = id_h, id_H_ = id_H_,
        ind_RE = object$model_data$ind_RE,
        W0_H = W0_H, W_H = W_H, U_H = U_H, X_H = X_H, Z_H = Z_H,
        Wlong_bar = object$Wlong_bar, Wlong_sds = object$Wlong_sds,
        any_gammas = any_gammas,
        FunForms_cpp = FunForms_cpp, FunForms_ind = FunForms_ind,
        Funs_FunForms = Funs_FunForms
    )

    H <- cum_haz(Data, components_newdata$mcmc)
    index <- rep(seq_along(times), n_times)
    for (i in seq_along(times)) {
        H[index == i, ] <- colCumsums(H[index == i, ])
    }
    CIF <- 1.0 - exp(- H)
    res <- list(pred = rowMeans(CIF),
                low = rowQuantiles(CIF, probs = (1 - level) / 2),
                upp = rowQuantiles(CIF, probs = (1 + level) / 2))
}


tt <- res#predict_Event(jointFit1, components_newdata, ND, level = 0.9)

matplot(times[[1]], cbind(res$low, res$pred, res$upp)[index == 2, ],
        type = "l", lty = c(2, 1, 2), col = c(1, 2, 1), lwd = 2, ylim = c(0, 1))



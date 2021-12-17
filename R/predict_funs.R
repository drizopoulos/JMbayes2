i_row <- function (x, i) x[i, ]

splt_REs <- function (b_mat, ind_RE) {
    n <- length(ind_RE)
    out <- vector("list", n)
    for (i in seq_len(n)) out[[i]] <- b_mat[, ind_RE[[i]], drop = FALSE]
    out
}

linpred_long <- function (X, betas, Z, b, id, type) {
    out <- vector("list", length(X))
    for (i in seq_along(X)) {
        out[[i]] <- c(X[[i]] %*% betas[[i]])
        if (type == "subject_specific") {
            out[[i]] <- out[[i]] +
                as.vector(rowSums(Z[[i]] * b[[i]][id[[i]], , drop = FALSE]))
        }
    }
    out
}

mu_fun <- function (eta, link) {
    switch (link, "identity" = eta, "inverse" = 1 / eta, "logit" = plogis(eta),
            "probit" = pnorm(eta), "cloglog" = - exp(- exp(eta)) + 1.0,
            "log" = exp(eta))
}

fix_NAs_preds <- function (preds, NAs, n) {
    if (is.null(NAs)) preds else {
        r <- rep(as.numeric(NA), n)
        r[-NAs] <- preds
        r
    }
}

prepare_Data_preds <- function (object, newdataL, newdataE) {
    # control
    control <- object$control
    # extract idVar and time_var
    idVar <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    # set dataL as newdata; almost the same code as in jm()
    # check for tibbles
    if (inherits(newdataL, "tbl_df") || inherits(newdataL, "tbl")) {
        newdataL <- as.data.frame(newdataL)
    }

    dataL <- newdataL
    idL <- dataL[[idVar]]
    nY <- length(unique(idL))
    # order data by idL and time_var
    if (is.null(dataL[[time_var]])) {
        stop("the variable specified in agument 'time_var' cannot be found ",
             "in the database of the longitudinal models.")
    }
    dataL <- dataL[order(idL, dataL[[time_var]]), ]

    # extract terms
    respVars <- object$model_info$var_names$respVars
    terms_FE <- object$model_info$terms$terms_FE
    terms_FE_noResp <- object$model_info$terms$terms_FE_noResp
    terms_RE <- object$model_info$terms$terms_RE
    terms_Surv <- object$model_info$terms$terms_Surv_noResp
    Xbar <- object$model_data$Xbar
    # create model frames
    mf_FE_dataL <- lapply(terms_FE, model.frame.default, data = dataL)
    mf_RE_dataL <- lapply(terms_RE, model.frame.default, data = dataL)

    # we need to account for missing data in the fixed and random effects model frames,
    # in parallel across outcomes (i.e., we will allow that some subjects may have no data
    # for some outcomes)
    NAs_FE_dataL <- lapply(mf_FE_dataL, attr, "na.action")
    NAs_RE_dataL <- lapply(mf_RE_dataL, attr, "na.action")
    mf_FE_dataL <- mapply2(fix_NAs_fixed, mf_FE_dataL, NAs_FE_dataL, NAs_RE_dataL)
    mf_RE_dataL <- mapply2(fix_NAs_random, mf_RE_dataL, NAs_RE_dataL, NAs_FE_dataL)

    # create response vectors
    y <- lapply(mf_FE_dataL, model.response)
    y <- lapply(y, function (yy) {
        if (is.factor(yy)) as.numeric(yy != levels(yy)[1L]) else yy
    })
    y[] <- lapply(y, as.matrix)
    NAs <- mapply2(c, NAs_FE_dataL, NAs_RE_dataL)
    times_y <- lapply(NAs, function (ind) if (!is.null(ind))
        dataL[[time_var]][-ind] else dataL[[time_var]])
    families <- object$model_info$families
    family_names <- sapply(families, "[[", "family")
    links <- sapply(families, "[[", "link")
    # for family = binomial and when y has two columns, set the second column
    # to the number of trials instead the number of failures
    binomial_data <- family_names %in% c("binomial", "beta binomial")
    trials_fun <- function (y) {
        if (NCOL(y) == 2L) y[, 2L] <- y[, 1L] + y[, 2L]
        y
    }
    y[binomial_data] <- lapply(y[binomial_data], trials_fun)
    unq_id <- unique(idL)
    idL <- mapply2(exclude_NAs, NAs_FE_dataL, NAs_RE_dataL,
                   MoreArgs = list(id = idL))
    idL <- lapply(idL, match, table = unq_id)
    idL_lp <- lapply(idL, function (x) match(x, unique(x)))
    unq_idL <- lapply(idL, unique)
    X <- mapply2(model.matrix.default, terms_FE, mf_FE_dataL)
    Z <- mapply2(model.matrix.default, terms_RE, mf_RE_dataL)

    ################################

    # extract terms
    terms_Surv <- object$model_info$terms$terms_Surv
    terms_Surv_noResp <- object$model_info$terms$terms_Surv_noResp
    type_censoring <- object$model_info$type_censoring
    # check for tibbles
    if (inherits(newdataE, "tbl_df") || inherits(newdataE, "tbl")) {
        newdataE <- as.data.frame(newdataE)
    }
    dataS <- newdataE
    idT <- dataS[[idVar]]
    mf_surv_dataS <- model.frame.default(terms_Surv, data = dataS)
    if (!is.null(NAs_surv <- attr(mf_surv_dataS, "na.action"))) {
        idT <- idT[-NAs_surv]
        dataS <- dataS[-NAs_surv, ]
    }
    idT <- factor(idT, levels = unique(idT))
    nT <- length(unique(idT))
    if (nY != nT) {
        stop("the number of groups/subjects in the longitudinal and survival datasets ",
             "do not seem to match. A potential reason why this may be happening is ",
             "missing data in some covariates used in the individual models.")
    }
    Surv_Response <- model.response(mf_surv_dataS)
    if (type_censoring == "right") {
        Time_right <- unname(Surv_Response[, "time"])
        Time_left <- Time_start <- trunc_Time <- rep(0.0, nrow(dataS))
        delta <-  unname(Surv_Response[, "status"])
    } else if (type_censoring == "counting") {
        Time_start <- unname(Surv_Response[, "start"])
        Time_stop <- unname(Surv_Response[, "stop"])
        delta <-  unname(Surv_Response[, "status"])
        Time_right <- Time_stop
        trunc_Time <- Time_start # possible left truncation time
        Time_left <- rep(0.0, nrow(dataS))
    } else if (type_censoring == "interval") {
        Time1 <-  unname(Surv_Response[, "time1"])
        Time2 <-  unname(Surv_Response[, "time2"])
        trunc_Time <- Time_start <- rep(0.0, nrow(dataS))
        delta <- unname(Surv_Response[, "status"])
        Time_right <- Time1
        Time_right[delta == 3] <- Time2[delta == 3]
        Time_right[delta == 2] <- 0.0
        Time_left <- Time1
        Time_left[delta <= 1] <- 0.0
    }
    if (type_censoring != "counting") {
        names(Time_right) <- names(Time_left) <- names(Time_start) <- idT
    }
    which_event <- which(delta == 1)
    which_right <- which(delta == 0)
    which_left <- which(delta == 2)
    which_interval <- which(delta == 3)
    # extract strata if present otherwise all subjects in one stratum
    ind_strata <- attr(terms_Surv, "specials")$strata
    strata <- if (is.null(ind_strata)) {
        rep(1, nrow(mf_surv_dataS))
    } else {
        unclass(mf_surv_dataS[[ind_strata]])
    }
    Time_integration <- Time_right
    Time_integration[which_left] <- Time_left[which_left]
    Time_integration[which_interval] <- Time_left[which_interval]
    Time_integration2 <- rep(0.0, length(Time_integration))
    if (length(which_interval)) {
        Time_integration2[which_interval] <- Time_right[which_interval]
    }
    last_times <- switch(type_censoring,
                         "right" = unname(Surv_Response[, "time"]),
                         "counting" = unname(Surv_Response[, "stop"]),
                         "interval" = unname(Surv_Response[, "time1"]))
    # create Gauss Kronrod points and weights
    GK <- gaussKronrod(control$GK_k)
    sk <- GK$sk
    P <- c(Time_integration - trunc_Time) / 2
    st <- outer(P, sk) + (c(Time_integration + trunc_Time) / 2)
    log_Pwk <- unname(rep(log(P), each = length(sk)) +
                          rep_len(log(GK$wk), length.out = length(st)))
    if (length(which_interval)) {
        # we take the absolute value because for the subjects for whom we do not have
        # interval censoring P2 will be negative and this will produce a NA when we take
        # the log in 'log_Pwk2'
        P2 <- abs(Time_integration2 - Time_integration) / 2
        st2 <- outer(P2, sk) + (c(Time_integration2 + Time_integration) / 2)
        log_Pwk2 <- rep(log(P2), each = length(sk)) +
            rep_len(log(GK$wk), length.out = length(st2))
    } else {
        P2 <- st2 <- log_Pwk2 <- rep(0.0, nT * control$GK_k)
    }

    # knots for the log baseline hazard function
    knots <- control$knots

    # indices
    ni_event <- tapply(idT, idT, length)
    ni_event <- cbind(c(0, head(cumsum(ni_event), -1)), cumsum(ni_event))
    id_H <- rep(paste0(idT, "_", unlist(tapply(idT, idT, seq_along))),
                each = control$GK_k)
    id_H <- match(id_H, unique(id_H))
    # id_H_ repeats each unique idT the number of quadrature points
    id_H_ <- rep(idT, each = control$GK_k)
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

    # Design matrices
    strata_H <- rep(strata, each = control$GK_k)
    W0_H <- create_W0(c(t(st)), knots, control$Bsplines_degree + 1, strata_H)
    idT_str <- paste0(idT, "_", strata)
    dataS_H <- SurvData_HazardModel(split(st, row(st)), dataS, Time_start,
                                    idT_str, time_var,
                                    match(idT_str, unique(idT_str)))
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
    if (length(which_event)) {
        W0_h <- create_W0(Time_right, knots, control$Bsplines_degree + 1,
                          strata)
        dataS_h <- SurvData_HazardModel(split(Time_right, seq_along(Time_right)),
                                        dataS, Time_start,
                                        idT_str, time_var,
                                        match(idT_str, unique(idT_str)))
        mf <- model.frame.default(terms_Surv_noResp, data = dataS_h)
        W_h <- construct_Wmat(terms_Surv_noResp, mf)
        if (!any_gammas) {
            W_h <- matrix(0.0, nrow = nrow(W_h), ncol = 1L)
        }
        X_h <- design_matrices_functional_forms(Time_right, terms_FE_noResp,
                                                dataL, time_var, idVar, idT,
                                                collapsed_functional_forms, Xbar,
                                                eps, direction)
        Z_h <- design_matrices_functional_forms(Time_right, terms_RE,
                                                dataL, time_var, idVar, idT,
                                                collapsed_functional_forms, NULL,
                                                eps, direction)
        U_h <- lapply(functional_forms, construct_Umat, dataS = dataS_h)
    } else {
        W0_h <- W_h <- matrix(0.0)
        X_h <- Z_h <- U_h <- rep(list(matrix(0.0)), length(respVars))
    }
    if (length(which_interval)) {
        W0_H2 <- create_W0(c(t(st2)), knots, control$Bsplines_degree + 1,
                           strata_H)
        dataS_H2 <- SurvData_HazardModel(split(st2, row(st2)), dataS, Time_start,
                                         idT_str, time_var,
                                         match(idT_str, unique(idT_str)))
        mf2 <- model.frame.default(terms_Surv_noResp, data = dataS_H2)
        W_h <- construct_Wmat(terms_Surv_noResp, mf2)
        if (!any_gammas) {
            W_H2 <- matrix(0.0, nrow = nrow(W_H2), ncol = 1L)
        }
        X_H2 <- design_matrices_functional_forms(st, terms_FE_noResp,
                                                 dataL, time_var, idVar, idT,
                                                 collapsed_functional_forms, Xbar,
                                                 eps, direction)
        Z_H2 <- design_matrices_functional_forms(st, terms_RE,
                                                 dataL, time_var, idVar, idT,
                                                 collapsed_functional_forms, NULL,
                                                 eps, direction)
        U_H2 <- lapply(functional_forms, construct_Umat, dataS = dataS_H2)
    } else {
        W0_H2 <- W_H2 <- matrix(0.0)
        X_H2 <- Z_H2 <- U_H2 <- rep(list(matrix(0.0)), length(respVars))
    }
    X_H[] <- lapply(X_H, docall_cbind)
    X_h[] <- lapply(X_h, docall_cbind)
    X_H2[] <- lapply(X_H2, docall_cbind)
    Z_H[] <- lapply(Z_H, docall_cbind)
    Z_h[] <- lapply(Z_h, docall_cbind)
    Z_H2[] <- lapply(Z_H2, docall_cbind)

    #W_bar <- object$W_bar
    #W_sds <- object$W_sds
    #W_H <- center_fun(W_H, W_bar, W_sds)
    #W_h <- center_fun(W_h, W_bar, W_sds)
    #W_H2 <- center_fun(W_H2, W_bar, W_sds)
    list(
        nY = nY, times_y = times_y, respVars = respVars, strata = strata,
        NAs_FE_dataL = NAs_FE_dataL, NAs_RE_dataL = NAs_RE_dataL,
        ind_RE = object$model_data$ind_RE,
        W0_H = W0_H, W0_h = W0_h, W0_H2 = W0_H2,
        W_H = W_H, W_h = W_h, W_H2 = W_H2,
        X_H = X_H, X_h = X_h, X_H2 = X_H2,
        Z_H = Z_H, Z_h = Z_h, Z_H2 = Z_H2,
        U_H = U_H, U_h = U_h, U_H2 = U_H2,
        Wlong_bar = object$Wlong_bar, Wlong_sds = object$Wlong_sds,
        idT = match(idT, unique(idT)), log_Pwk = log_Pwk, log_Pwk2 = log_Pwk2,
        id_H = id_H, id_H_ = id_H_, id_h = id_h, any_gammas = any_gammas,
        which_event = which_event, which_right = which_right,
        which_left = which_left, which_interval = which_interval,
        ni_event = ni_event, FunForms_cpp = FunForms_cpp,
        FunForms_ind = FunForms_ind, Funs_FunForms = Funs_FunForms,
        X = X, Z = Z, y = y, family_names = family_names,
        links = links, extra_parms = object$model_data$extra_parms,
        unq_idL = unq_idL, idL_lp = idL_lp, idL = idL, last_times = last_times,
        Time_start = Time_start
    )
}

prepare_DataE_preds <- function (object, newdataL, newdataE,
                                 low_limit, upp_limit, last_times, n_times,
                                 st0 = NULL, index = NULL, index2 = NULL) {
    # control
    control <- object$control
    # extract idVar and time_var
    idVar <- object$model_info$var_names$idVar
    respVars <- unlist(object$model_info$var_names$respVars)
    time_var <- object$model_info$var_names$time_var
    # set dataL as newdata; almost the same code as in jm()
    # check for tibbles
    if (inherits(newdataL, "tbl_df") || inherits(newdataL, "tbl")) {
        newdataL <- as.data.frame(newdataL)
    }

    dataL <- newdataL
    idL <- dataL[[idVar]]
    # order data by idL and time_var
    if (is.null(dataL[[time_var]])) {
        stop("the variable specified in agument 'time_var' cannot be found ",
             "in the database of the longitudinal models.")
    }
    dataL <- dataL[order(idL, dataL[[time_var]]), ]
    terms_FE_noResp <- object$model_info$terms$terms_FE_noResp
    terms_RE <- object$model_info$terms$terms_RE
    Xbar <- object$model_data$Xbar

    ################################

    # extract terms
    terms_Surv <- object$model_info$terms$terms_Surv
    terms_Surv_noResp <- object$model_info$terms$terms_Surv_noResp
    # check for tibbles
    if (inherits(newdataE, "tbl_df") || inherits(newdataE, "tbl")) {
        newdataE <- as.data.frame(newdataE)
    }
    dataS <- newdataE
    idT <- dataS[[idVar]]
    mf_surv_dataS <- model.frame.default(terms_Surv, data = dataS)
    if (!is.null(NAs_surv <- attr(mf_surv_dataS, "na.action"))) {
        idT <- idT[-NAs_surv]
        dataS <- dataS[-NAs_surv, ]
    }
    idT <- factor(idT, levels = unique(idT))
    nT <- length(unique(idT))
    Surv_Response <- model.response(mf_surv_dataS)
    delta <- unname(Surv_Response[, "status"])
    which_event <- which(delta == 1)
    which_right <- which(delta == 0)
    which_left <- which(delta == 2)
    which_interval <- which(delta == 3)
    # extract strata if present otherwise all subjects in one stratum
    ind_strata <- attr(terms_Surv, "specials")$strata
    strata <- if (is.null(ind_strata)) {
        rep(1, nrow(mf_surv_dataS))
    } else {
        unclass(mf_surv_dataS[[ind_strata]])
    }
    # create Gauss Kronrod points and weights
    GK <- gaussKronrod(control$GK_k)
    sk <- GK$sk
    P <- c(upp_limit - low_limit) / 2
    st <- outer(P, sk) + (c(upp_limit + low_limit) / 2)
    log_Pwk <- unname(rep(log(P), each = length(sk)) +
                          rep_len(log(GK$wk), length.out = length(st)))
    P2 <- st2 <- log_Pwk2 <- rep(0.0, nT * control$GK_k)

    # knots for the log baseline hazard function
    knots <- control$knots

    # indices
    ni_event <- tapply(idT, idT, length)
    ni_event <- cbind(c(0, head(cumsum(ni_event), -1)), cumsum(ni_event))
    id_H <- c(t(row(st)))
    id_H <- match(id_H, unique(id_H))
    # id_H_ repeats each unique idT the number of quadrature points
    id_H_ <- rep(idT, each = control$GK_k)
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

    # Design matrices
    strata_H <- rep(strata, each = control$GK_k)
    W0_H <- create_W0(c(t(st)), knots, control$Bsplines_degree + 1, strata_H)
    dataS_H <- SurvData_HazardModel(split(st, row(st)), dataS, last_times,
                                    paste0(idT, "_", strata), time_var,
                                    rep(index, each = control$GK_k))
    mf <- model.frame.default(terms_Surv_noResp, data = dataS_H)
    W_H <- construct_Wmat(terms_Surv_noResp, mf)
    any_gammas <- as.logical(ncol(W_H))
    if (!any_gammas) {
        W_H <- matrix(0.0, nrow = nrow(W_H), ncol = 1L)
    }
    attr <- lapply(functional_forms, extract_attributes, data = dataS_H)
    eps <- lapply(attr, "[[", 1L)
    direction <- lapply(attr, "[[", 2L)
    if (is.null(index2)) {
        index2 <- match(idT, unique(idT))
    }
    if (object$model_info$CR_MS) {
        index2_H <- rep(index2, each = control$GK_k)
    } else {
        index2_H <- rep(index2, n_times)
    }
    X_H <- design_matrices_functional_forms(split(st, row(st)), terms_FE_noResp,
                                            dataL, time_var, idVar, index2_H,
                                            collapsed_functional_forms, Xbar,
                                            eps, direction)
    Z_H <- design_matrices_functional_forms(split(st, row(st)), terms_RE,
                                            dataL, time_var, idVar,index2_H,
                                            collapsed_functional_forms, NULL,
                                            eps, direction)
    U_H <- lapply(functional_forms, construct_Umat, dataS = dataS_H)
    if (length(which_event)) {
        W0_h <- create_W0(c(t(st0)), knots, control$Bsplines_degree + 1,
                          strata)
        dataS_h <- SurvData_HazardModel(split(st0, row(st0)), dataS, last_times,
                                        paste0(idT, "_", strata), time_var,
                                        index)
        mf <- model.frame.default(terms_Surv_noResp, data = dataS_h)
        W_h <- construct_Wmat(terms_Surv_noResp, mf)
        if (!any_gammas) {
            W_h <- matrix(0.0, nrow = nrow(W_h), ncol = 1L)
        }
        X_h <- design_matrices_functional_forms(split(st0, row(st0)), terms_FE_noResp,
                                                dataL, time_var, idVar, index2,
                                                collapsed_functional_forms, Xbar,
                                                eps, direction)
        Z_h <- design_matrices_functional_forms(split(st0, row(st0)), terms_RE,
                                                dataL, time_var, idVar, index2,
                                                collapsed_functional_forms, NULL,
                                                eps, direction)
        U_h <- lapply(functional_forms, construct_Umat, dataS = dataS_h)
    } else {
        W0_h <- W_h <- matrix(0.0)
        X_h <- Z_h <- U_h <- rep(list(matrix(0.0)), length(respVars))
    }
    if (length(which_interval)) {
        W0_H2 <- create_W0(c(t(st2)), knots, control$Bsplines_degree + 1,
                           strata_H)
        dataS_H2 <- SurvData_HazardModel(split(st2, row(st2)), dataS, last_times,
                                         paste0(idT, "_", strata), time_var,
                                         rep(index, each = control$GK_k))
        mf2 <- model.frame.default(terms_Surv_noResp, data = dataS_H2)
        W_h <- construct_Wmat(terms_Surv_noResp, mf2)
        if (!any_gammas) {
            W_H2 <- matrix(0.0, nrow = nrow(W_H2), ncol = 1L)
        }
        X_H2 <- design_matrices_functional_forms(split(st2, row(st2)), terms_FE_noResp,
                                                 dataL, time_var, idVar,
                                                 rep(index2, each = control$GK_k),
                                                 collapsed_functional_forms, Xbar,
                                                 eps, direction)
        Z_H2 <- design_matrices_functional_forms(split(st2, row(st2)), terms_RE,
                                                 dataL, time_var, idVar,
                                                 rep(index2, each = control$GK_k),
                                                 collapsed_functional_forms, NULL,
                                                 eps, direction)
        U_H2 <- lapply(functional_forms, construct_Umat, dataS = dataS_H2)
    } else {
        W0_H2 <- W_H2 <- matrix(0.0)
        X_H2 <- Z_H2 <- U_H2 <- rep(list(matrix(0.0)), length(respVars))
    }
    X_H[] <- lapply(X_H, docall_cbind)
    X_h[] <- lapply(X_h, docall_cbind)
    X_H2[] <- lapply(X_H2, docall_cbind)
    Z_H[] <- lapply(Z_H, docall_cbind)
    Z_h[] <- lapply(Z_h, docall_cbind)
    Z_H2[] <- lapply(Z_H2, docall_cbind)
    list(
        ind_RE = object$model_data$ind_RE,
        W0_H = W0_H, W0_h = W0_h, W0_H2 = W0_H2,
        W_H = W_H, W_h = W_h, W_H2 = W_H2,
        X_H = X_H, X_h = X_h, X_H2 = X_H2,
        Z_H = Z_H, Z_h = Z_h, Z_H2 = Z_H2,
        U_H = U_H, U_h = U_h, U_H2 = U_H2,
        Wlong_bar = object$Wlong_bar, Wlong_sds = object$Wlong_sds,
        idT = match(idT, unique(idT)), log_Pwk = log_Pwk, log_Pwk2 = log_Pwk2,
        id_H = id_H, id_H_ = id_H_, id_h = id_h, any_gammas = any_gammas,
        which_event = which_event, which_right = which_right,
        which_left = which_left, which_interval = which_interval,
        ni_event = ni_event, FunForms_cpp = FunForms_cpp,
        FunForms_ind = FunForms_ind, Funs_FunForms = Funs_FunForms,
        strata = strata)
}

get_components_newdata <- function (object, newdata, n_samples, n_mcmc,
                                    cores, seed) {
    if (!exists(".Random.seed", envir = .GlobalEnv)) {
        runif(1L)
    }
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))

    # prepare the data for calculations
    newdataL <- if (!is.data.frame(newdata)) newdata[["newdataL"]] else newdata
    newdataE <- if (!is.data.frame(newdata)) newdata[["newdataE"]] else newdata
    if (!object$model_info$CR_MS) {
        idT <- newdataE[[idVar <- object$model_info$var_names$idVar]]
        keep_last_row <- tapply(row.names(newdataE), factor(idT, unique(idT)), tail, 1L)
        newdataE <- newdataE[keep_last_row, ]
    }
    Data <- prepare_Data_preds(object, newdataL, newdataE)

    # MCMC sample
    b <- lapply(sapply(Data$Z, ncol), function (nc) matrix(0.0, Data$nY, nc))
    M <- sum(sapply(object$mcmc$bs_gammas, nrow))
    get_param <- function (nam) {
        tht <- object$mcmc[[nam]]
        if (!is.null(tht)) docall_rbind(tht) else matrix(0.0, M, 1)
    }
    ind_betas <- grep("^betas", names(object$mcmc))
    mcmc <- list(
        b = b, bs_gammas = get_param("bs_gammas"), gammas = get_param("gammas"),
        alphas = get_param("alphas"),
        Wlong_std_alphas = get_param("Wlong_std_alphas"),
        W_std_gammas = get_param("W_std_gammas"),
        betas = lapply(object$mcmc[ind_betas], docall_rbind)
    )
    has_sigmas <- object$model_data$has_sigmas
    mcmc$sigmas <- matrix(0.0, M, length(has_sigmas))
    mcmc$sigmas[, has_sigmas > 0] <- get_param("sigmas")
    D <- get_param("D")
    mcmc$D <- array(0.0, c(dim(lowertri2mat(D[1L, ])), M))
    for (i in seq_len(M)) {
        mcmc$D[, , i] <- lowertri2mat(D[i, ])
    }
    if (n_samples > M) {
        warning("the number of samples cannot be greater than the number of ",
                "MCMC iterations in the fitted model.")
        n_samples <- M
    }
    control <- list(GK_k = object$control$GK_k, n_samples = n_samples, n_iter = n_mcmc)
    id_samples <-
        split(seq_len(control$n_samples),
              rep(seq_len(cores), each = ceiling(control$n_samples / cores),
                  length.out = control$n_samples))

    sample_parallel <- function (id_samples, Data, mcmc, control) {
        # keep only the samples from the MCMC used in the sampling
        # of the random effects
        mcmc$bs_gammas <- mcmc$bs_gammas[id_samples, , drop = FALSE]
        mcmc$gammas <- mcmc$gammas[id_samples, , drop = FALSE]
        mcmc$alphas <- mcmc$alphas[id_samples, , drop = FALSE]
        mcmc$betas[] <- lapply(mcmc$betas, function (m, ind) m[ind, , drop = FALSE],
                               "ind" = id_samples)
        mcmc$sigmas <- mcmc$sigmas[id_samples, , drop = FALSE]
        mcmc$D <- mcmc$D[, , id_samples, drop = FALSE]
        # update control n_samples
        control$n_samples <- length(id_samples)
        # update random effects
        mcmc[["b"]] <- simulate_REs(Data, mcmc, control)
        mcmc$Wlong_std_alphas <- mcmc$Wlong_std_alphas[id_samples, , drop = FALSE]
        mcmc$W_std_gammas <- mcmc$W_std_gammas[id_samples, , drop = FALSE]
        mcmc
    }
    if (cores > 1L) {
        cl <- parallel::makeCluster(cores)
        parallel::clusterSetRNGStream(cl = cl, iseed = seed)
        out <- parallel::parLapply(cl, id_samples, sample_parallel,
                                   Data = Data, mcmc = mcmc, control = control)
        parallel::stopCluster(cl)
    } else {
        set.seed(seed)
        out <- list(sample_parallel(id_samples[[1L]], Data = Data, mcmc = mcmc,
                                    control = control))
    }
    combine <- function (x) {
        n <- length(x)
        res <- x[[1L]]
        if (n > 1L) {
            for (i in 2:n) {
                res$bs_gammas <- rbind(res$bs_gammas, x[[i]][["bs_gammas"]])
                res$gammas <- rbind(res$gammas, x[[i]][["gammas"]])
                res$alphas <- rbind(res$alphas, x[[i]][["alphas"]])
                res$sigmas <- rbind(res$sigmas, x[[i]][["sigmas"]])
                res$Wlong_std_alphas <-
                    rbind(res$Wlong_std_alphas, x[[i]][["Wlong_std_alphas"]])
                res$W_std_gammas <-
                    rbind(res$W_std_gammas, x[[i]][["W_std_gammas"]])
                d1 <- dim(res$D)[3L]
                d2 <- dim(x[[i]][["D"]])[3L]
                a <- array(0.0, dim = c(dim(res$D)[1:2], d1 + d2))
                a[, , seq(1, d1)] <- res$D
                a[, , seq(d1 + 1, d1 + d2)] <- x[[i]][["D"]]
                res$D <- a
                d1 <- dim(res$b)[3L]
                d2 <- dim(x[[i]][["b"]])[3L]
                a <- array(0.0, dim = c(dim(res$b)[1:2], d1 + d2))
                a[, , seq(1, d1)] <- res$b
                a[, , seq(d1 + 1, d1 + d2)] <- x[[i]][["b"]]
                res$b <- a
                for (j in seq_along(res$betas)) {
                    res$betas[[j]] <- rbind(res$betas[[j]], x[[i]][["betas"]][[j]])
                }
            }
        }
        res
    }
    list(mcmc = combine(out), X = Data$X, Z = Data$Z, y = Data$y,
         times_y = Data$times_y, id = Data$idL,
         ind_RE = Data$ind_RE, links = Data$links,
         respVars = lapply(Data$respVars, "[", 1L),
         NAs = mapply2(c, Data$NAs_FE_dataL, Data$NAs_RE_dataL),
         last_times = Data$last_times, strata = Data$strata, idT = Data$idT)
}

predict_Long <- function (object, components_newdata, newdata, newdata2, times,
                          type, type_pred, level, return_newdata, return_mcmc) {
    # Predictions for newdata
    newdataL <- if (!is.data.frame(newdata)) newdata[["newdataL"]] else newdata
    betas <- components_newdata$mcmc[["betas"]]
    b_mat <- components_newdata$mcmc[["b"]]
    ind_RE <- components_newdata$ind_RE
    links <- components_newdata$links
    K <- length(ind_RE)
    M <- dim(b_mat)[3L]
    out <- lapply(components_newdata$X, function (x) matrix(0.0, nrow(x), M))
    names(out) <- components_newdata$respVars
    for (i in seq_len(M)) {
        eta_i <-
            linpred_long(components_newdata$X, lapply(betas, i_row, i),
                         components_newdata$Z,
                         splt_REs(rbind(b_mat[, , i]), ind_RE),
                         components_newdata$id, type = type)
        for (j in seq_len(K)) {
            out[[j]][, i] <- if (type_pred == "response") {
                mu_fun(eta_i[[j]], links[j])
            } else eta_i[[j]]
        }
    }
    res1 <- list(preds = lapply(out, rowMeans, na.rm = TRUE),
                 low = lapply(out, rowQuantiles, probs = (1 - level) / 2),
                 upp = lapply(out, rowQuantiles, probs = (1 + level) / 2),
                 mcmc = if (return_mcmc) out)
    if (return_newdata) {
        n <- nrow(newdataL)
        preds <- mapply2(fix_NAs_preds, res1$preds, components_newdata$NAs,
                         MoreArgs = list(n = n))
        names(preds) <- paste0("pred_", components_newdata$respVars)
        low <- mapply2(fix_NAs_preds, res1$low, components_newdata$NAs,
                       MoreArgs = list(n = n))
        names(low) <- paste0("low_", components_newdata$respVars)
        upp <- mapply2(fix_NAs_preds, res1$upp, components_newdata$NAs,
                       MoreArgs = list(n = n))
        names(upp) <- paste0("upp_", components_newdata$respVars)
        l <- c(preds, low, upp)
        l <- l[c(matrix(seq_along(l), ncol = length(preds), byrow = TRUE))]
        res1 <- cbind(newdataL, as.data.frame(do.call("cbind", l)))
    }
    ############################################################################
    ############################################################################
    # Predictions for newdata2
    if (is.null(newdata2) && !is.null(times) && is.numeric(times)) {
        last_times <- components_newdata$last_times
        if (object$model_info$CR_MS) {
            ff <- components_newdata$idT
            ff <- factor(ff, levels = unique(ff))
            last_times <- tapply(last_times, ff, max)
        }
        t_max <- max(object$model_data$Time_right)
        test <- sapply(last_times, function (lt, tt) all(tt <= lt), tt = times)
        if (any(test)) {
            stop("according to the definition of argument 'times', for some ",
                 "subjects the last available time is\n\t larger than the ",
                 "maximum time to predict; redefine 'times' accordingly.")
        }
        f <- function (lt, tt, tm) c(lt, sort(tt[tt > lt & tt <= tm]))
        times <- lapply(last_times, f, tt = times, tm = t_max)
        n_times <- sapply(times, length)
        newdata2 <- newdataL
        idVar <- object$model_info$var_names$idVar
        time_var <- object$model_info$var_names$time_var
        idT <- newdata2[[idVar]]
        newdata2 <- newdata2[tapply(row.names(newdata2),
                                    factor(idT, unique(idT)), tail, 1L), ]
        newdata2 <- newdata2[rep(seq_along(times), n_times), ]
        newdata2[[time_var]] <- unlist(times, use.names = FALSE)
    }
    if (!is.null(newdata2)) {
        newdataL2 <- if (!is.data.frame(newdata2)) newdata2[["newdataL"]] else newdata2
        terms_FE_noResp <- object$model_info$terms$terms_FE_noResp
        terms_RE <- object$model_info$terms$terms_RE
        mf_FE <- lapply(terms_FE_noResp, model.frame.default, data = newdataL2)
        mf_RE <- lapply(terms_RE, model.frame.default, data = newdataL2)
        NAs_FE <- lapply(mf_FE, attr, "na.action")
        NAs_RE <- lapply(mf_RE, attr, "na.action")
        mf_FE <- mapply2(fix_NAs_fixed, mf_FE, NAs_FE, NAs_RE)
        mf_RE <- mapply2(fix_NAs_random, mf_RE, NAs_RE, NAs_FE)
        X <- mapply2(model.matrix.default, terms_FE_noResp, mf_FE)
        Z <- mapply2(model.matrix.default, terms_RE, mf_RE)
        NAs <- mapply2(c, NAs_FE, NAs_RE)
        idL <- newdataL2[[object$model_info$var_names$idVar]]
        unq_id <- unique(idL)
        idL <- mapply2(exclude_NAs, NAs_FE, NAs_RE, MoreArgs = list(id = idL))
        idL <- lapply(idL, match, table = unq_id)
        out <- lapply(X, function (x) matrix(0.0, nrow(x), M))
        names(out) <- components_newdata$respVars
        for (i in seq_len(M)) {
            eta_i <-
                linpred_long(X, lapply(betas, i_row, i), Z,
                             splt_REs(rbind(b_mat[, , i]), ind_RE),
                             idL, type = type)
            for (j in seq_len(K)) {
                out[[j]][, i] <- if (type_pred == "response") {
                    mu_fun(eta_i[[j]], links[j])
                } else eta_i[[j]]
            }
        }
        res2 <- list(preds = lapply(out, rowMeans, na.rm = TRUE),
                     low = lapply(out, rowQuantiles, probs = (1 - level) / 2),
                     upp = lapply(out, rowQuantiles, probs = (1 + level) / 2),
                     mcmc = if (return_mcmc) out)
        if (return_newdata) {
            n <- nrow(newdataL2)
            preds <- mapply2(fix_NAs_preds, res2$preds, NAs,
                             MoreArgs = list(n = n))
            names(preds) <- paste0("pred_", components_newdata$respVars)
            low <- mapply2(fix_NAs_preds, res2$low, NAs,
                           MoreArgs = list(n = n))
            names(low) <- paste0("low_", components_newdata$respVars)
            upp <- mapply2(fix_NAs_preds, res2$upp, NAs,
                           MoreArgs = list(n = n))
            names(upp) <- paste0("upp_", components_newdata$respVars)
            l <- c(preds, low, upp)
            l <- l[c(matrix(seq_along(l), ncol = length(preds), byrow = TRUE))]
            res2 <- cbind(newdataL2, as.data.frame(do.call("cbind", l)))
        }
    }
    out <- if (is.null(newdata2)) {
        res1
    } else {
        list(newdata = res1, newdata2 = res2)
    }
    class(out) <- c("predict_jm", class(out))
    attr(out, "id_var") <- object$model_info$var_names$idVar
    attr(out, "time_var") <- object$model_info$var_names$time_var
    attr(out, "resp_vars") <- object$model_info$var_names$respVars_form
    attr(out, "ranges") <- ranges <- lapply(object$model_data$y, range,
                                            na.rm = TRUE)
    attr(out, "last_times") <-
        with(components_newdata, tapply(last_times, idT, tail, 1L))
    attr(out, "y") <- components_newdata$y
    attr(out, "times_y") <- components_newdata$times_y
    attr(out, "id") <- components_newdata$id
    attr(out, "process") <- "longitudinal"
    out
}

predict_Event <- function (object, components_newdata, newdata, newdata2,
                           times, level, return_newdata, return_mcmc) {
    # prepare the data for calculations
    newdataL <- if (!is.data.frame(newdata)) newdata[["newdataL"]] else newdata
    newdataE <- if (!is.data.frame(newdata)) newdata[["newdataE"]] else newdata
    if (!is.null(newdata2)) {
        newdataE2 <- if (!is.data.frame(newdata2)) newdata2[["newdataE"]] else newdata2
    } else {
        newdataE2 <- newdataE
    }
    CR_MS <- object$model_info$CR_MS

    # The definition of last_times needs to be checked for counting and interval
    last_times <- components_newdata$last_times
    if (CR_MS) {
        last_times2 <- last_times
        ff <- paste(components_newdata$idT, components_newdata$strata, sep = "_")
        ff <- factor(ff, levels = unique(ff))
        last_times <- tapply(last_times, ff, tail, n = 1L)
    }
    t_max <- quantile(object$model_data$Time_right, probs = 0.9)
    if (is.null(times) || !is.numeric(times)) {
        g <- function (tt, t_max) {
            if (tt > t_max) seq(tt, max(object$model_data$Time_right), length.out = 21L)
            else seq(tt, t_max, length.out = 21L)
        }
        times <- lapply(last_times, g, t_max = t_max)
    } else {
        t_max <- max(object$model_data$Time_right)
        test <- sapply(last_times, function (lt, tt) all(tt <= lt), tt = times)
        if (any(test)) {
            stop("according to the definition of argument 'times', for some ",
                 "subjects the last available time is \nlarger than the ",
                 "maximum time to predict; redefine 'times' accordingly.")
        }
        f <- function (lt, tt, tm) c(lt, sort(tt[tt > lt & tt <= tm]))
        times <- lapply(last_times, f, tt = times, tm = t_max)
    }
    n_times <- sapply(times, length)

    if (!CR_MS) {
        idT <- newdataE[[idVar <- object$model_info$var_names$idVar]]
        idT <- factor(idT, unique(idT))
        keep_last_row <- tapply(row.names(newdataE), idT, tail, 1L)
        newdataE <- newdataE[keep_last_row, ]
        upp_limit <- unlist(times, use.names = FALSE)
        gg <- function (t0, t) c(t0, head(t, -1))
        low_limit <- unlist(mapply2(gg, last_times, times), use.names = FALSE)
        Data <- prepare_DataE_preds(object, newdataL, newdataE, low_limit,
                                    upp_limit, last_times, n_times)
        n_strata <- length(unique(Data$strata))
        idt_str <- rep(unique(Data$idT), each = n_strata)
        n_times_id <- tapply(n_times, idt_str, sum)
        Data$id_H_ <- rep(unique(idt_str), n_times_id * object$control$GK_k)
        H <- cum_haz(Data, components_newdata$mcmc)
        index <- rep(seq_along(times), n_times)
        for (i in seq_along(times)) {
            H[index == i, ] <- colCumsums(H[index == i, ])
        }
        CIF <- 1.0 - pmax(exp(- H), .Machine$double.eps)
        newdataE <- newdataE[rep(seq_along(times), n_times), ]
        newdataE[[object$model_info$var_names$time_var]] <-
            unlist(times, use.names = FALSE)
    } else {
        Data <- prepare_Data_preds(object, newdataL, newdataE)
        idT <- Data$idT
        H <- cum_haz(Data, components_newdata$mcmc)
        S0 <- exp(- rowsum(H, idT, reorder = FALSE))
        upp_limit <- unlist(times, use.names = FALSE)
        gg <- function (t0, t) c(t0, head(t, -1))
        low_limit <- unlist(mapply2(gg, last_times, times), use.names = FALSE)
        GK <- gaussKronrod(object$control$GK_k)
        sk <- GK$sk
        P <- c(upp_limit - low_limit) / 2
        st <- outer(P, sk) + (c(upp_limit + low_limit) / 2)
        log_Pwk <- unname(rep(log(P), each = length(sk)) +
                              rep_len(log(GK$wk), length.out = length(st)))
        n_strata <- length(unique(Data$strata))
        upp_limit2 <- c(t(st))
        low_limit2 <- 0 * upp_limit2
        last_row_str <- tapply(row.names(newdataE2), ff, tail, n = 1L)
        newdataE2[[object$model_info$var_names$event_var]] <- 1
        Data2 <- prepare_DataE_preds(object, newdataL, newdataE2,
                                     low_limit2, upp_limit2, last_times2,
                                     n_times, st0 = st,
                                     index = rep(seq_along(n_times), n_times),
                                     index2 = rep(rep(seq_along(unique(idT)),
                                                      each = n_strata), n_times))
        idt_str <- rep(unique(Data2$idT), each = n_strata)
        n_times_id <- tapply(n_times, idt_str, sum)
        Data2$id_H <- rep(seq_along(st), each = object$control$GK_k)
        Data2$id_H_ <- rep(unique(idt_str), n_times_id * object$control$GK_k^2)
        Data2$id_h <- rep(idt_str, n_times * object$control$GK_k)
        Data2$which_event <- seq_len(nrow(Data2$W0_h))

        log_hS <- hSfun(Data2, components_newdata$mcmc)
        ind <- c(t(row(st)))
        hS <- rowsum.default(exp(log_Pwk + log_hS), ind, reorder = FALSE)

        index <- rep(seq_along(times), n_times)
        for (i in seq_along(times)) {
            hS[index == i, ] <- colCumsums(hS[index == i, ])
        }
        CIF <- hS / S0[rep(unique(idt_str), n_times_id), ]
        newdataE2 <- newdataE2[last_row_str, ]
        newdataE2 <- newdataE2[rep(seq_along(times), n_times), ]
        newdataE2[[object$model_info$var_names$time_var]] <-
            unlist(times, use.names = FALSE)
    }
    res <- list(pred = rowMeans(CIF),
                low = rowQuantiles(CIF, probs = (1 - level) / 2),
                upp = rowQuantiles(CIF, probs = (1 + level) / 2),
                times = unlist(times, use.names = FALSE),
                id = rep(levels(idT), n_times),
                "_strata" = if (CR_MS)
                    rep(tapply(components_newdata$strata, ff, tail, 1), n_times),
                mcmc = if (return_mcmc) CIF
            )
    if (return_newdata) {
        newdataE2[["pred_CIF"]] <- res$pred
        newdataE2[["low_CIF"]] <- res$low
        newdataE2[["upp_CIF"]] <- res$upp
        newdataE2[["_strata"]] <- res[["_strata"]]
        res <- newdataE2
    }
    class(res) <- c("predict_jm", class(res))
    attr(res, "id_var") <- object$model_info$var_names$idVar
    attr(res, "time_var") <- object$model_info$var_names$time_var
    attr(res, "resp_vars") <- object$model_info$var_names$respVars_form
    attr(res, "ranges") <- ranges <- lapply(object$model_data$y, range,
                                            na.rm = TRUE)
    attr(res, "last_times") <-
        with(components_newdata, tapply(last_times, idT, tail, 1L))
    attr(res, "y") <- components_newdata$y
    attr(res, "times_y") <- components_newdata$times_y
    attr(res, "id") <- components_newdata$id
    attr(res, "process") <- "event"
    res
}




jm <- function (Surv_object, Mixed_objects, time_var, recurrent = FALSE,
                functional_forms = NULL, which_independent = NULL,
                base_hazard = NULL, data_Surv = NULL, id_var = NULL,
                priors = NULL, control = NULL, ...) {
    call <- match.call()
    # control argument:
    # - GK_k: number of quadrature points for the Gauss Kronrod rule; options 15 and 7
    # - Bsplines_degree: the degree of the splines in each basis; default quadratic splines
    # - base_hazard_segments: number of segments to split the follow-up period; default 10
    # - diff: the order of the difference used in the penalty matrix for the
    #         B-splines for h_0; default is 2
    # - n_chains: the number of chains for the MCMC
    # - n_burnin: the number of burn-in iterations
    # - n_iter: the number of total iterations per chain
    # - seed: the seed used in the sampling procedures
    # - parallel: what type of parallel computing to use, "snow" (default) or
    #             "multicore"
    # - cores: the number of cores to use for running the chains in parallel;
    #          no point of setting this greater than n_chains
    # - knots: the knots for the log baseline hazard B-spline approximation
    # - MALA: if TRUE, the MALA algorithm is used to update the elements of
    #         of the Cholesky factor of the D matrix
    # - save_random_effects: if TRUE, the random effects are stored in the fitted object
    # - knots: numeric vector with the knots for the baseline hazard function
    con <- list(GK_k = 15L, n_chains = 3L, n_burnin = 500L, n_iter = 3500L,
                n_thin = 1L, seed = 123L, MALA = FALSE,
                save_random_effects = FALSE, save_logLik_contributions = FALSE,
                basis = "bs", Bsplines_degree = 2L, base_hazard_segments = 9L,
                knots = NULL, timescale_base_hazard = "identity", diff = 2L,
                parallel = "snow",
                cores = parallelly::availableCores(omit = 1L))
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0) {
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    }
    if (con$n_burnin > con$n_iter) {
        stop("'n_burnin' cannot be larger than 'n_iter'.")
    }
    # if a single mixed model has been provided put in a list
    if (!inherits(Mixed_objects, "list")) {
        Mixed_objects <- list(Mixed_objects)
    }
    # check if only lme and MixMod have been provided
    if (!all(sapply(Mixed_objects, class) %in% c("lme", "MixMod"))) {
        stop("'Mixed_objects' should be of class 'lme' of 'MixMod'.\n")
    }
    # extract the data from each of the mixed models
    # and check whether the same data have been used;
    # otherwise an error
    datas <- lapply(Mixed_objects, "[[", "data")
    datas[] <- lapply(datas, function (d)
        if (inherits(d, "tbl_df") || inherits(d, "tbl")) as.data.frame(d) else d)
    if (!all(sapply(datas[-1L], function (x) isTRUE(all.equal(x, datas[[1L]]))))) {
        stop("It seems that some of the mixed models have been fitted to different versions ",
             "of the dataset. Use the same exact dataset in the calls to lme() ",
             " and mixed_model().")
    }
    dataL <- datas[[1L]]
    rm(datas)
    # extract id variable (again we assume a single grouping variable)
    id_names <- sapply(Mixed_objects, function (object)
        names(if (inherits(object, "MixMod")) object$id[1L] else object$groups[1L]))
    if (!all(id_names == id_names[1L])) {
        stop("it seems that different grouping variables have been used in the mixed models.")
    }
    idVar <- id_names[1L]
    idL <- dataL[[idVar]]
    idL <- factor(idL, levels = unique(idL))
    nY <- length(unique(idL))
    # order data by idL and time_var
    if (is.null(dataL[[time_var]])) {
        stop("the variable specified in agument 'time_var' cannot be found ",
             "in the database of the longitudinal models.")
    }
    dataL <- dataL[order(idL, dataL[[time_var]]), ]
    idL <- dataL[[idVar]]

    # extract terms from mixed models
    terms_FE <- lapply(Mixed_objects, extract_terms, which = "fixed", data = dataL)
    respVars <- lapply(terms_FE, function (tt) all.vars(attr(tt, "variables")[[2L]]))
    respVars_form <- sapply(terms_FE, function (tt) as.character(attr(tt, "variables"))[2L])
    terms_FE_noResp <- lapply(terms_FE, delete.response)
    terms_RE <- lapply(Mixed_objects, extract_terms, which = "random", data = dataL)

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
    if (any(sapply(y, function (x) any(!is.finite(x))))) {
        stop("infite value detected in some longitudinal outcomes. These are not allowed.\n")
    }

    # extract families
    families <- lapply(Mixed_objects, "[[", "family")
    families[sapply(families, is.null)] <- rep(list(gaussian()),
                                               sum(sapply(families, is.null)))

    # extra parameter in the families
    extra_parms <- sapply(families,
                          function (x) if (is.null(x$df)) 0.0 else x$df)

    # create the idL per outcome
    # IMPORTANT: some ids may be missing when some subjects have no data for a particular outcome
    # This needs to be taken into account when using idL for indexing. Namely, a new id variable
    # will need to be created in jm_fit()
    unq_id <- unique(idL)
    idL <- mapply2(exclude_NAs, NAs_FE_dataL, NAs_RE_dataL,
                   MoreArgs = list(id = idL))
    idL[] <- lapply(idL, match, table = unq_id)
    # the index variable idL_lp is to be used to subset the random effects of each outcome
    # such that to calculate the Zb part of the model as rowSums(Z * b[idL_lp, ]). This
    # means that for outcomes that miss some subjects, we recode the id variable from 1
    # until n', where n' is the number of subjects available for the respective outcome
    idL_lp <- lapply(idL, function (x) match(x, unique(x)))
    # the unique values of idL is to be used in specifying which subjects have which outcomes
    # this is relevant in the calculation of the log density / probability mass function
    # for the longitudinal outcomes
    unq_idL <- lapply(idL, unique)

    # create design matrices for mixed models
    X <- mapply2(model.matrix.default, terms_FE, mf_FE_dataL)
    Xbar <- lapply(X, colMeans)
    Z <- mapply2(model.matrix.default, terms_RE, mf_RE_dataL)
    if (length(Z) == 1 && ncol(Z[[1]]) == 1) {
        stop("jm() does not currently work when you have a single ",
             "longitudinal outcome and only random intercepts.")
    }
    nres <- sapply(Z, ncol)
    ind_RE <- split(seq_len(sum(nres)), rep(seq_along(Z), nres))
    componentsHC <- mapply2(create_HC_X, x = X, z = Z, id = idL,
                            terms = terms_FE, data = mf_FE_dataL)
    x_in_z <- lapply(componentsHC, "[[", "x_in_z")
    x_notin_z <- lapply(componentsHC, "[[", "x_notin_z")
    x_in_z_base <- lapply(componentsHC, "[[", "x_in_z_base")
    xbas_in_z <- lapply(componentsHC, "[[", "xbas_in_z")
    z_in_x <- lapply(componentsHC, "[[", "z_in_x")
    Xbar[] <- mapply2(center_X, Xbar, x_notin_z)
    X_HC <- lapply(componentsHC, "[[", "X_HC")
    mat_HC <- lapply(componentsHC, "[[", "mat_HC")
    nfes <- sapply(X, ncol)
    # 'ind_FE' is used in vec2field() to re-create the field of betas
    # from betas_vec
    ind_FE <- split(seq_len(sum(nfes)), rep(seq_along(X), nfes))
    # 'ind_FE_HC' denotes which elements of betas_vec are in the HC formulation
    # this will be used to save the results in the corresponding columns
    ind_FE_HC <- unlist(mapply2(function (x, ind) x[ind], ind_FE, x_in_z_base),
                        use.names = FALSE)
    #
    # 'ind_FE_HC' denotes which elements of betas_vec are not in the
    # HC formulation. It is a list, to be used in conjunction with
    # has_tilde_betas. That is, if has_tilde_betas = TRUE, we need to save
    # from the Metropolis-Hastings step the betas in the columns 'ind_FE_nHC'
    ind_FE_nHC <- mapply2(function (x, ind) x[-ind], ind_FE, x_in_z_base)
    has_tilde_betas <- as.integer(sapply(ind_FE_nHC, length) > 0)
    ind_FE_nHC[] <- lapply(ind_FE_nHC, function (x) if (length(x)) x else 0L)
    if (any(!unlist(lapply(x_notin_z, is.na), use.names = FALSE)) && !any(namc %in% "n_iter")) {
        con$n_iter <- 6500L
    }
    ########################################################
    ########################################################
    # try to recover survival dataset
    if (is.null(data_Surv) || !is.data.frame(data_Surv)) {
        dataS <- try(eval(Surv_object$call$data, envir = parent.frame()),
                     silent = TRUE)
        if (inherits(dataS, "try-error")) {
            stop("could not recover the dataset used to fit the Cox/AFT model; ",
                 "please provide this dataset in the 'data_Surv' argument of jm().")
        }
    } else {
        dataS <- data_Surv
    }
    if (inherits(dataS, "tbl_df") || inherits(dataS, "tbl")) {
        dataS <- as.data.frame(dataS)
    }
    # if the longitudinal outcomes are not in dataS, we set a random value for
    # them. This is needed for the calculation of the matrix of interaction terms
    # between the longitudinal outcomes and other variables.
    for (i in seq_along(respVars)) {
        nl <- length(respVars[[i]])
        for (j in seq_len(nl)) {
            if (is.null(dataS[[respVars[[i]][j]]]))
                dataS[[respVars[[i]][j]]] <- rnorm(nrow(dataS), 20)
        }
    }
    # if the time_var is not in dataS set it to a random number
    if (is.null(dataS[[time_var]])) dataS[[time_var]] <- rnorm(nrow(dataS))
    # terms for survival model
    terms_Surv <- Surv_object$terms
    terms_Surv_noResp <- delete.response(terms_Surv)
    mf_surv_dataS <- model.frame.default(terms_Surv, data = dataS)
    # var names
    av <- all.vars(attr(terms_Surv, "variables")[[2L]])
    Time_var <- head(av, -1L)
    event_var <- tail(av, 1L)
    # survival times
    if (!is.null(NAs_surv <- attr(mf_surv_dataS, "na.action"))) {
        dataS <- dataS[-NAs_surv, ]
    }
    Surv_Response <- model.response(mf_surv_dataS)
    type_censoring <- attr(Surv_Response, "type")
    if (is.null(dataS[[idVar]])) {
        if (is.null(id_var)) {
            stop("cannot extract the subject id variable from the dataset used to fit the ",
                 "survival model. Please include this variable in the dataset ",
                 "and/or specify the 'id_var' argument.\n")
        } else {
            idT <- dataS[[id_var]]
        }
    } else {
        idT <- dataS[[idVar]]
    }
    idT <- factor(idT, levels = unique(idT))

    nT <- length(unique(idT))
    if (nY != nT) {
        stop("the number of groups/subjects in the longitudinal and survival datasets ",
             "do not seem to match. A potential reason why this may be happening is ",
             "missing data in some covariates used in the individual models.")
    }
    if (!all(idT %in% dataL[[idVar]])) {
        stop("it seems that some of the levels of the id variable in the survival dataset",
             "cannot be found in the dataset of the Mixed_objects. Please check that ",
             "the same subjects/groups are used in the datasets used to fit the mixed ",
             "and survival models. Also, the name of the subjects/groups variable ",
             "in the different datasets used to fit the individual models ",
             "needs to be the same in all of the datasets.")
    }
    # we need to check that the ordering of the subjects in the same in dataL and dataS.
    # If not, then a warning and do it internally
    if (!all(order(unique(idT)) == order(factor(unq_id, levels = unq_id)))) {
        warning("It seems that the ordering of the subjects in the dataset used to fit the ",
                "mixed models and the dataset used for the survival model is not the same. ",
                "We set internally the datasets in the same order, but it would be best ",
                "that you do it beforehand on your own.")
        dataS <- dataS[order(idT), ]
        idT <- dataS[[id_var]]
        idT <- factor(idT, levels = unique(idT))
        mf_surv_dataS <- model.frame.default(terms_Surv, data = dataS)
        Surv_Response <- model.response(mf_surv_dataS)
    }
    # Notation:
    #  - Time_right: event or right censoring time
    #  - Time_left: left censoring time
    #  - trunc_Time: truncation time
    #  - delta: 0 of right censored, 1 for event, 2 for left censored,
    #           3 for interval censored
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
        if(!"GK_k" %in% namc) con$GK_k <- 7L
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
        strt <- mf_surv_dataS[[ind_strata]]
        if (!is.factor(strt)) {
            stop("the strata variable in the Cox model is not a factor. ",
                 "Please set it as a factor in the database,\n\tand ",
                 "refit the stratified Cox model.")
        }
        unclass(strt)
    }
    n_strata <- length(unique(strata))
    con$basis <- rep_len(con$basis, n_strata)
    con$Bsplines_degree <- rep_len(con$Bsplines_degree, n_strata)
    con$base_hazard_segments <- rep_len(con$base_hazard_segments, n_strata)
    con$timescale_base_hazard <- rep_len(con$timescale_base_hazard, n_strata)
    con$diff <- rep_len(con$diff, n_strata)
    if (!is.null(base_hazard)) {
        if (length(base_hazard) != n_strata) {
            base_hazard <- rep_len(base_hazard, n_strata)
            warning("length of 'base_hazard' has been replicated to match the ",
                    "length of strata.\n")
        }
        bhh <- matrix(NA_real_, n_strata, 6)
        for (k in seq_len(n_strata)) {
            if (!is.na(base_hazard[k])) {
                bhh[k, ] <- bh <- get_hazard(base_hazard[k])
                if (bh['log_time']) con$timescale_base_hazard[k] <- "log"
                if (bh['ns']) con$basis[k] <- "ns"
                if (bh['pwc_const']) con$Bsplines_degree[k] <- 0L
                if (bh['pwc_linear']) con$Bsplines_degree[k] <- 1L
                if (bh['weibull']) {
                    con$timescale_base_hazard[k] <- "log"
                    con$basis[k] <- "ns"
                    con$base_hazard_segments[k] <- 1L
                }
            }
        }
    }

    # extract weights if present, otherwise set weight equal to one for all subjects
    indx <- match("weights", names(Surv_object$call), nomatch = 0)
    weights <- if (indx) {
        if (is.null(Surv_object$model)) {
            stop("Please refit the model using coxph(..., model = TRUE).\n")
        }
        model.weights(Surv_object$model)
    } else {
        rep(1, nrow(mf_surv_dataS))
    }
    intgr_ind <- attr(weights, "integrate")
    if (is.null(intgr_ind)) intgr_ind <- 0
    intgr <- any(as.logical(intgr_ind))

    # check if we have competing risks or multi-state processes. In this case,
    # we will have multiple strata per subject. NOTE: this will also be the
    # case for recurrent events. We will need to change the definition of
    # CR_MS to account for this
    CR_MS <- any(tapply(strata, idT, function (x) length(unique(x))) > 1)

    # 'Time_integration' is the upper limit of the integral in likelihood
    # of the survival model. For subjects with event (delta = 1), for subjects with
    # right censoring and for subjects with interval censoring we need to integrate
    # up to 'Time_right'. For subjects with left censoring we need to integrate up to
    # 'Time_left'; hence we set for them 'Time_integration = Time_left'.
    # For subjects with interval censoring we need two integrals from 0 to 'Time_right'
    # and also from 0 to 'Time_left'. However, in the Gauss-Kronrod approximation it
    # can happen that the first integral has a lower value than the second one, which is
    # not correct. To overcome this issue, for interval censored data we first approximate
    # the integral from 0 to 'Time_left'. And then the integral from 0 to 'Time_right' is
    # set equal to first integral plus the integral from 'Time_left' to 'Time_right'. For
    # this second integral we introduce the variable 'Time_integration2' which is equal to
    # 'Time_right', i.e., for interval censored data 'Time_integration' is set equal to
    # 'Time_left'
    Time_integration <- Time_right
    Time_integration[which_left] <- Time_left[which_left]
    Time_integration[which_interval] <- Time_left[which_interval]
    Time_integration2 <- rep(0.0, length(Time_integration))
    if (length(which_interval)) {
        Time_integration2[which_interval] <- Time_right[which_interval]
    }

    # create Gauss Kronrod points and weights
    GK <- gaussKronrod(con$GK_k)
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
        P2 <- st2 <- log_Pwk2 <- rep(0.0, nT * con$GK_k)
    }

    # knots for the log baseline hazard function
    if (is.null(con$knots)) {
        qs <- if (recurrent == "gap") {
            lapply(split(Time_right - trunc_Time, strata), range)
        } else {
            lapply(split(Time_right, strata), range)
        }
        con$knots <- mapply2(knots, x = qs, ndx = con$base_hazard_segments,
                             deg = con$Bsplines_degree, basis = con$basis)
    }
    for (k in seq_len(n_strata)) {
        if (con$timescale_base_hazard[k] != "identity") {
            con$knots[[k]] <- log(con$knots[[k]])
        }
    }

    # Extract functional forms per longitudinal outcome
    if (is.language(functional_forms)) {
        term_labels <- attr(terms(functional_forms), "term.labels")
        ind_tlabs <- lapply(respVars_form, grep, term_labels, fixed = TRUE)
        functional_forms <- lapply(ind_tlabs, function (ind, tlabs)
            if (length(ind)) reformulate(tlabs[ind]), tlabs = term_labels)
        names(functional_forms) <- respVars_form
        functional_forms <- functional_forms[!sapply(functional_forms, is.null)]
    }
    if (any(!names(functional_forms) %in% respVars_form)) {
        stop("unknown names in the list provided in the 'functional_forms' argument; as names ",
             "of the elements of this list you need to use the response variables from ",
             "the multivariate mixed model.\n")
    }
    # for outcomes not specified in Formulas use the value parameterization
    not_specified <- !respVars_form %in% names(functional_forms)
    if (any(not_specified)) {
        functional_forms_ns <- lapply(respVars_form[not_specified],
                                      function (v) as.formula(paste0("~ value(", v, ")")))
        names(functional_forms_ns) <- respVars_form[not_specified]
        functional_forms <- c(functional_forms, functional_forms_ns)
    }
    functional_forms <- functional_forms[order(match(names(functional_forms),
                                                     respVars_form))]
    functional_forms <- mapply2(expand_Dexps, functional_forms, respVars_form)

    ###################################################################
    # List of lists
    # One list component per association structure per outcome
    # List components vectors of integers corresponding to the term
    # each association structure corresponds to
    FunForms_per_outcome <-
        mapply2(extract_functional_forms, Form = functional_forms,
                nam = respVars_form, MoreArgs = list(data = dataS))
    FunForms_per_outcome <- lapply(FunForms_per_outcome,
                                   function (x) x[sapply(x, length) > 0])

    collapsed_functional_forms <- lapply(FunForms_per_outcome, names)

    collapsed_functional_forms <-
        lapply(collapsed_functional_forms, function (nam) {
            nn <- c("value", "slope", "area", "velocity", "acceleration",
                    "coefs", "Delta")
            names(unlist(sapply(nn, grep, x = nam, fixed = TRUE,
                                simplify = FALSE)))
        })

    Funs_FunForms <-
        mapply2(extractFuns_FunForms, Form = functional_forms, nam = respVars_form,
                MoreArgs = list(data = dataS))
    #####################################################

    # design matrices for the survival submodel:
    #  - W0 is the design matrix for the log baseline hazard
    #  - W is the design matrix for the covariates in the Surv_object
    #    (including exogenous time-varying covariates)
    #  - X is the design matrix for the fixed effects, per outcome and functional form
    #  - Z is the design matrix for the random effects, per outcome and functional form
    #  - U is the design matrix for possible interaction terms in functional forms
    #  - Wlong is the design matrix for the longitudinal outcomes in the survival submodel
    #    that is already multiplied with the interaction terms matrix U
    # in the above design matrices we put the "_h" to denote calculation at the event time
    # 'Time_right', we put "_H" to denote calculation at the 'Time_integration', and
    # "_H2" to denote calculation at the 'Time_integration2'.
    strata_H <- rep(strata, each = con$GK_k)
    W0_H <- if (recurrent == "gap") {
        create_W0(c(t(st - trunc_Time)), con$knots, con$Bsplines_degree,
                  strata_H, con$basis, con$timescale_base_hazard)
    } else {
        create_W0(c(t(st)), con$knots, con$Bsplines_degree, strata_H,
                  con$basis, con$timescale_base_hazard)
    }
    dataS_H <- SurvData_HazardModel(split(st, row(st)), dataS, Time_start,
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
    zero_ind <- lapply(attr, "[[", 3L)
    zero_ind_X <- lapply(zero_ind, function (x) if (length(x)) x[[1]][["X"]] else x)
    zero_ind_Z <- lapply(zero_ind, function (x) if (length(x)) x[[1]][["Z"]] else x)
    time_window <- lapply(attr, "[[", 4L)
    standardise <- lapply(attr, "[[", 5L)
    IE_time <- IE_time2 <- lapply(attr, "[[", 6L) #!! new
    if(any(lengths(IE_time))) { #!! new
      bool <- vapply(IE_time, function(x)
        all(vapply(x, function(x_i) identical(x_i, x[[1]]), logical(1))), 
        logical(1))
      if (any(!bool)) {
        bad <- names(IE_time)[!bool]
        msg <- "Inconsistent 'IE_time' variable detected within in the following outcome(s): "
        stop(paste0(msg, paste(bad, collapse = ", ")), ".")
      }
      remove(bool)
      IE_time <- lapply(IE_time, function(x) if (length(x)) x[[1L]] else x)
      if (!all(vapply(IE_time, function(x) is.character(x) || !length(x), logical(1)))) {
        stop("'IE_time' must be a character string.")
      }
      missing_IE_time <- setdiff(unlist(IE_time), colnames(dataS))
      if (length(missing_IE_time)) {
        stop("Cannot extract the IE_time variable(s) from the dataset used to ",
             "fit the survival model. Please include the following variable(s)",
             " in the dataset: ", paste(missing_IE_time, collapse = ", "), ".")
      }
      remove(missing_IE_time)
    }
    IE_time <- lapply(IE_time, function(v) if (length(v)) dataS[[v]] else v) #!! new
    X_H <- design_matrices_functional_forms(st, terms_FE_noResp,
                                            dataL, time_var, idVar, idT,
                                            collapsed_functional_forms, Xbar,
                                            eps, direction, zero_ind_X, 
                                            time_window, standardise, IE_time) #!! new
    Z_H <- design_matrices_functional_forms(st, terms_RE,
                                            dataL, time_var, idVar, idT,
                                            collapsed_functional_forms, NULL,
                                            eps, direction, zero_ind_Z, 
                                            time_window, standardise, IE_time) #!! new
    U_H <- lapply(functional_forms, construct_Umat, dataS = dataS_H)
    if (length(which_event)) {
        W0_h <- if (recurrent == "gap") {
            create_W0(Time_right - trunc_Time, con$knots,
                      con$Bsplines_degree, strata, con$basis,
                      con$timescale_base_hazard)
        } else {
            create_W0(Time_right, con$knots, con$Bsplines_degree,
                      strata, con$basis, con$timescale_base_hazard)
        }
        dataS_h <- SurvData_HazardModel(split(Time_right, seq_along(Time_right)),
                                        dataS, Time_start,
                                        paste0(idT, "_", strata), time_var)
        mf <- model.frame.default(terms_Surv_noResp, data = dataS_h)
        W_h <- construct_Wmat(terms_Surv_noResp, mf)
        if (!any_gammas) {
            W_h <- matrix(0.0, nrow = nrow(W_h), ncol = 1L)
        }
        X_h <- design_matrices_functional_forms(Time_right, terms_FE_noResp,
                                                dataL, time_var, idVar, idT,
                                                collapsed_functional_forms, Xbar,
                                                eps, direction, zero_ind_X, 
                                                time_window, standardise, IE_time) #!! new
        Z_h <- design_matrices_functional_forms(Time_right, terms_RE,
                                                dataL, time_var, idVar, idT,
                                                collapsed_functional_forms, NULL,
                                                eps, direction, zero_ind_Z, 
                                                time_window, standardise, IE_time) #!! new
        U_h <- lapply(functional_forms, construct_Umat, dataS = dataS_h)
    } else {
        W0_h <- W_h <- matrix(0.0)
        X_h <- Z_h <- U_h <- rep(list(matrix(0.0)), length(respVars))
    }
    if (length(which_interval)) {
        W0_H2 <- create_W0(c(t(st2)), con$knots, con$Bsplines_degree, strata_H,
                      con$basis, con$timescale_base_hazard)
        dataS_H2 <- SurvData_HazardModel(split(st2, row(st2)), dataS, Time_start,
                                         paste0(idT, "_", strata), time_var)
        mf2 <- model.frame.default(terms_Surv_noResp, data = dataS_H2)
        W_H2 <- construct_Wmat(terms_Surv_noResp, mf2)
        if (!any_gammas) {
            W_H2 <- matrix(0.0, nrow = nrow(W_H2), ncol = 1L)
        }
        X_H2 <- design_matrices_functional_forms(st2, terms_FE_noResp,
                                                 dataL, time_var, idVar, idT,
                                                 collapsed_functional_forms, Xbar,
                                                 eps, direction, zero_ind_X, 
                                                 time_window, standardise, 
                                                 IE_time) #!! new
        Z_H2 <- design_matrices_functional_forms(st2, terms_RE,
                                                 dataL, time_var, idVar, idT,
                                                 collapsed_functional_forms, NULL,
                                                 eps, direction, zero_ind_Z, 
                                                 time_window, standardise, 
                                                 IE_time) #!! new
        U_H2 <- lapply(functional_forms, construct_Umat, dataS = dataS_H2)
    } else {
        W0_H2 <- W_H2 <- matrix(0.0)
        X_H2 <- Z_H2 <- U_H2 <- rep(list(matrix(0.0)), length(respVars))
    }
    nfes_HC <- sapply(x_in_z_base, length)
    out_in <- sapply(idL, "%in%", x = seq_len(nT))
    all_pat <- apply(out_in, 1L, paste0, collapse = "/")
    id_patt <- match(all_pat, unique(all_pat))
    find_patt <- function (patt, n) which(rep(patt, times = n))
    unq_out_in <- split(unique(out_in), seq_len(nrow(unique(out_in))))
    ind_RE_patt <- lapply(unq_out_in, find_patt, n = nres)
    ind_FE_patt <- lapply(unq_out_in, find_patt, n = nfes_HC)
    X_dot <- create_X_dot(nres, nfes_HC, z_in_x, x_in_z, X_HC, nT, unq_idL,
                          xbas_in_z)
    ############################################################################
    ############################################################################
    Data <- list(n = nY, idL = idL, idL_lp = idL_lp, unq_idL = unq_idL,
                 y = y, X = X, Z = Z, X_dot = X_dot, Xbar = Xbar,
                 x_in_z = x_in_z, x_notin_z = x_notin_z,
                 has_tilde_betas = has_tilde_betas, ind_FE = ind_FE,
                 ind_FE_HC = ind_FE_HC, ind_FE_nHC = ind_FE_nHC,
                 id_patt = id_patt, ind_RE_patt = ind_RE_patt,
                 ind_FE_patt = ind_FE_patt,
                 #####
                 idT = idT, any_gammas = any_gammas, strata = strata,
                 weights = weights, intgr = intgr, intgr_ind = as.numeric(intgr_ind),
                 Time_right = Time_right, Time_left = Time_left, Time_start = Time_start,
                 delta = delta, which_event = which_event, which_right = which_right,
                 which_left = which_left, which_interval = which_interval,
                 W0_H = W0_H, W_H = W_H, X_H = X_H, Z_H = Z_H, U_H = U_H,
                 W0_h = W0_h, W_h = W_h, X_h = X_h, Z_h = Z_h, U_h = U_h,
                 W0_H2 = W0_H2, W_H2 = W_H2, X_H2 = X_H2, Z_H2 = Z_H2, U_H2 = U_H2,
                 log_Pwk = log_Pwk, log_Pwk2 = log_Pwk2,
                 ind_RE = ind_RE, extra_parms = extra_parms,
                 which_term_h = lapply(seq_len(n_strata)[-1], function(x) which(x == strata)),
                 which_term_H = lapply(seq_len(n_strata)[-1], function(x) which(x == strata_H)))
    ############################################################################
    ############################################################################
    # objects to export
    data <- list(dataL = dataL, dataS = dataS)
    model_info <- list(
        terms = list(terms_FE = terms_FE, terms_FE_noResp = terms_FE_noResp,
                     terms_RE = terms_RE, terms_Surv = terms_Surv,
                     terms_Surv_noResp = terms_Surv_noResp),
        frames = list(mf_FE = mf_FE_dataL, mf_RE = mf_RE_dataL,
                      mf_Surv = mf_surv_dataS),
        var_names = list(respVars = respVars, respVars_form = respVars_form,
                         idVar = idVar, time_var = time_var,
                         Time_var = Time_var, event_var = event_var),
        families = families,
        type_censoring = type_censoring, CR_MS = CR_MS,
        functional_forms = functional_forms,
        FunForms_per_outcome = FunForms_per_outcome,
        collapsed_functional_forms = collapsed_functional_forms,
        FunForms_cpp = lapply(FunForms_per_outcome, unlist),
        FunForms_ind = FunForms_ind(FunForms_per_outcome),
        Funs_FunForms = lapply(Funs_FunForms, function (x) if (!is.list(x)) list(x) else x),
        eps = eps, direction = direction, zero_ind_X = zero_ind_X,
        zero_ind_Z = zero_ind_Z, time_window = time_window,
        standardise = standardise, IE_time = IE_time2, #!! new
        recurrent = !isFALSE(recurrent),
        ind_RE_patt = ind_RE_patt, ind_FE_patt = ind_FE_patt,
        id_patt = id_patt, callS = Surv_object$call,
        NAs_FE = NAs_FE_dataL, NAs_RE = NAs_RE_dataL
    )
    ############################################################################
    ############################################################################
    # initial values
    betas <- lapply(Mixed_objects, fixef)
    log_sigmas <- sapply(Mixed_objects, extract_log_sigmas)
    Data$has_sigmas <- as.integer(log_sigmas > -20)
    # indicator denoting subject-specific sigmas
    ss_sigmas <- rep(FALSE, length(log_sigmas))
    Data$ss_sigmas <- ss_sigmas
    sigmas <- split(exp(log_sigmas), seq_along(log_sigmas))
    sigmas[ss_sigmas] <-
        mapply2(rep, x = sigmas[ss_sigmas],
                length.out = lapply(idL_lp, lng_unq)[ss_sigmas])
    D_lis <- lapply(Mixed_objects, extract_D)
    D <- bdiag2(D_lis, which_independent = which_independent)
    if (abs(max(nearPD(D) - D)) > sqrt(.Machine$double.eps)) {
        D <- bdiag2(D_lis, off_diag_val = 1e-04,
                    which_independent = which_independent)
        D <- nearPD(D)
    }
    ind_zero_D <- which(abs(D) < sqrt(.Machine$double.eps), arr.ind = TRUE)
    ind_zero_D <- ind_zero_D[ind_zero_D[, 'col'] > ind_zero_D[, 'row'], , drop = FALSE]
    Data$ind_zero_D <- ind_zero_D
    b <- mapply2(extract_b, Mixed_objects, unq_idL,
                 MoreArgs = list(n = nY, unq_id = as.character(unq_id)))
    bs_gammas <- rep(-0.1, ncol(W0_H))
    gammas <- if (inherits(Surv_object, "coxph")) coef(Surv_object) else
        -coef(Surv_object)[-1L] / Surv_object$scale
    if (is.null(gammas)) gammas <- 0.0
    alphas <- rep(0.0, sum(sapply(U_H, ncol)))
    frailty <- rep(0.0, nT)
    alphaF <- rep(0.0, max(n_strata - 1, 1))
    sigmaF <- 0.1
    initial_values <-
        list(betas = betas, log_sigmas = log_sigmas,
             sigmas = sigmas, D = D, b = b, bs_gammas = bs_gammas,
             gammas = gammas, alphas = alphas,
             tau_bs_gammas = ifelse(con$base_hazard_segments > 1L, 20, 0.001),
             alphaF = alphaF, frailty = frailty, sigmaF = sigmaF)
    ############################################################################
    ############################################################################
    # Limits
    b_mat <- docall_cbind(b)
    Data$limit_b_low <- b_mat - rep(25 * sqrt(diag(D)), each = nrow(b_mat))
    Data$limit_b_upp <- b_mat + rep(25 * sqrt(diag(D)), each = nrow(b_mat))
    ############################################################################
    ############################################################################
    # Priors
    Tau_bs_gammas <-
        mapply2(function (nc, d)
            crossprod(diff(diag(nc), differences = if (nc > 3) d else 1)),
                nc = attr(W0_H, "ncW0"), d = con$diff)
    mean_betas <- lapply(betas, unname)
    mean_betas <- mapply2(function (b, m) {b[1] <- b[1] + m; b},
                          mean_betas, mapply2("%*%", Xbar, betas))
    Tau_betas <- mapply2(weak_informative_Tau, Mixed_objects, Xbar)
    #mean_betas <- lapply(betas, "*", 0.0)
    #Tau_betas <- lapply(betas, function (b) 0.01 * diag(length(b)))
    mean_betas_HC <- unlist(mean_betas, use.names = FALSE)[ind_FE_HC]
    Tau_betas_HC <- bdiag(Tau_betas)[ind_FE_HC, ind_FE_HC, drop = FALSE]
    mean_betas_nHC <- mapply2(get_betas_nHC, mean_betas, x_notin_z)
    Tau_betas_nHC <- mapply2(get_betas_nHC, Tau_betas, x_notin_z)
    prs <- list(mean_betas_HC = mean_betas_HC, Tau_betas_HC = Tau_betas_HC,
                mean_betas_nHC = mean_betas_nHC, Tau_betas_nHC = Tau_betas_nHC,
                mean_bs_gammas = lapply(Tau_bs_gammas, function (x) x[, 1] * 0),
                penalized_bs_gammas = rep(TRUE, n_strata),
                Tau_bs_gammas = Tau_bs_gammas,
                A_tau_bs_gammas = rep(5, n_strata),
                B_tau_bs_gammas = rep(0.5, n_strata),
                rank_Tau_bs_gammas =
                    sapply(lapply(Tau_bs_gammas, qr), "[[", 'rank'),
                mean_gammas = gammas,
                Tau_gammas = diag(0.25, length(gammas)),
                penalty_gammas = "none",
                A_lambda_gammas = 0.5, B_lambda_gammas = 1,
                A_tau_gammas = 0.5, B_tau_gammas = 1,
                A_nu_gammas = 0.5, B_nu_gammas = 1,
                A_xi_gammas = 0.5, B_xi_gammas = 1,
                mean_alphas = lapply(alphas, "*", 0.0),
                Tau_alphas = lapply(alphas, function (a) diag(0.25, length(a))),
                penalty_alphas = "none",
                A_lambda_alphas = 0.5, B_lambda_alphas = 1.0,
                A_tau_alphas = 0.5, B_tau_alphas = 1.0,
                A_nu_alphas = 0.5, B_nu_alphas = 1.0,
                A_xi_alphas = 0.5, B_xi_alphas = 1.0,
                gamma_prior_D_sds = TRUE,
                D_sds_df = 3.0, D_sds_sigma = rep(2.5, nrow(D)),
                D_sds_shape = 5.0, D_sds_mean = sqrt(diag(D)),
                D_L_etaLKJ = 3.0, gamma_prior_sigmas = TRUE,
                sigmas_df = 3.0,
                sigmas_sigmas = rep(5.0, length(log_sigmas)),
                sigmas_shape = 5.0, sigmas_mean = exp(log_sigmas),
                mean_alphaF = lapply(alphaF, "*", 0.0),
                Tau_alphaF = rep(list(diag(0.25, 1)), length(alphaF)),
                gamma_prior_sigmaF = TRUE,
                sigmaF_df = 3.0,
                sigmaF_sigmas = 5.0,
                sigmaF_shape = 0.25/0.4,
                sigmaF_mean = 0.25)
    if (!is.null(base_hazard)) {
        weib_ind <- !is.na(bhh[, 4]) & bhh[, 4] > 0
        prs$penalized_bs_gammas[weib_ind] <- FALSE
    }
    prs$Tau_bs_gammas[!prs$penalized_bs_gammas] <-
        lapply(prs$Tau_bs_gammas[!prs$penalized_bs_gammas],
               function (m) 0.1 *  diag(, ncol(m)))
    if (is.null(priors) || !is.list(priors)) {
        priors <- prs
    } else {
        ind <- intersect(names(priors), names(prs))
        prs[ind] <- priors[ind]
        priors <- prs
    }
    if (!priors$penalty_gammas %in% c("none", "ridge", "horseshoe")) {
        warning("'priors$penalty_gammas' can only take values, 'none', ",
                "'single' or 'double'.")
        priors$penalty_gammas <- "none"
    }
    if (priors$penalty_gammas != "none" && length(initial_values$gammas) == 1) {
        warning("it is more meaningful to penalize the 'gamma' coefficients ",
                "if their length is more than one.")
        priors$penalty_gammas <- "none"
    }
    if (!priors$penalty_alphas %in% c("none", "ridge", "horseshoe")) {
        warning("'priors$penalty_alphas' can only take values, 'none', ",
                "'single' or 'double'.")
        priors$penalty_alphas <- "none"
    }
    if (priors$penalty_alphas != "none" &&
        length(unlist(initial_values$alphas)) == 1) {
        warning("it is more meaningful to penalize the 'alpha' coefficients ",
                "if their length is more than one.")
        priors$penalty_alphas <- "none"
    }
    if (priors$penalty_gammas == "ridge") {
        priors$A_lambda_gammas <- 1
        priors$B_lambda_gammas <- 1
    }
    initial_values$tau_bs_gammas[!priors$penalized_bs_gammas] <- 1
    ############################################################################
    ############################################################################
    # Fit the model
    Fit <- jm_fit(Data, model_info, initial_values, priors, con)
    out <- c(Fit, list(model_data = c(Data, data), model_info = model_info,
                       initial_values = initial_values, control = con,
                       priors = priors, call = call))
    class(out) <- "jm"
    out
}

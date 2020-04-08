jm <- function (Surv_object, Mixed_objects, time_var,
                functional_forms = NULL, data_Surv = NULL, id_var = NULL,
                priors = NULL, control = NULL, ...) {
    call <- match.call()
    # control argument:
    # - GK_k: number of quadrature points for the Gauss Kronrod rule; options 15 and 7
    # - Bsplines_degree: the degree of the splines in each basis; default quadratic splines
    # - base_hazard_segments: number of segments to split the follow-up period; default 10
    # - n_chains: the number of chains for the MCMC
    # - n_adapt: the number of iterations to use in optimizing acceptance rates;
    #            this will be also the burn-in
    # - n_iter: the number of iterations per chain. These will be the iterations after
    #           n_adapt
    con <- list(GK_k = 15L, Bsplines_degree = 2, base_hazard_segments = 10,
                n_chains = 3L, n_adapt = 500L, n_iter = 1000L)
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0)
        warning("unknown names in control: ", paste(noNms, collapse = ", "))

    # extract the data from each of the mixed models
    # and check whether the same data have been used;
    # otherwise an error
    datas <- lapply(Mixed_objects, "[[", "data")
    if (!all(sapply(datas[-1L], all.equal, datas[[1L]]))) {
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
    nY <- length(unique(idL))
    # order data by idL and time_var
    dataL <- dataL[order(idL, dataL[[time_var]]), ]

    # extract terms from mixed models
    terms_FE <- lapply(Mixed_objects, extract_terms, which = "fixed", data = dataL)
    respVars <- sapply(terms_FE, function (tt) all.vars(tt)[1L])
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
    mf_FE_dataL <- mapply(fix_NAs_fixed, mf_FE_dataL, NAs_FE_dataL, NAs_RE_dataL)
    mf_RE_dataL <- mapply(fix_NAs_random, mf_RE_dataL, NAs_RE_dataL, NAs_FE_dataL)

    # create response vectors
    y <- lapply(mf_FE_dataL, model.response)
    y <- lapply(y, function (yy) {
        if (is.factor(yy)) as.numeric(yy != levels(yy)[1L]) else yy
    })
    if (any(sapply(y, function (x) any(!is.finite(x))))) {
        stop("infite value detected in some longitudinal outcomes. These are not allowed.\n")
    }

    # exctract families
    families <- lapply(Mixed_objects, "[[", "family")
    families[sapply(families, is.null)] <- rep(list(gaussian()),
                                               sum(sapply(families, is.null)))
    # create the idL per outcome
    # IMPORTANT: some ids may be missing when some subjects have no data for a particular outcome
    # This needs to be taken into account when using idL for indexing. Namely, a new id variable
    # will need to be created in jm_fit()
    unq_id <- unique(idL)
    idL <- mapply(exclude_NAs, NAs_FE_dataL, NAs_RE_dataL,
                  MoreArgs = list(id = idL), SIMPLIFY = FALSE)
    idL <- lapply(idL, match, table = unq_id)
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
    X <- mapply(model.matrix.default, terms_FE, mf_FE_dataL)
    Z <- mapply(model.matrix.default, terms_RE, mf_RE_dataL)
    componentsHC <- mapply(create_HC_X, terms_FE, terms_RE, X, Z, idL, mf_FE_dataL,
                           SIMPLIFY = FALSE)
    Xhc <- lapply(componentsHC, "[[", "Xhc")
    columns_HC <- lapply(componentsHC, "[[", "columns_HC")
    columns_nHC <- lapply(componentsHC, "[[", "columns_nHC")
    ########################################################
    ########################################################
    # try to recover survival dataset
    if (is.null(data_Surv))
        try(dataS <- eval(Surv_object$call$data, envir = parent.frame()),
            silent = TRUE)
    if (inherits(dataS, "try-error")) {
        stop("could not recover the dataset used to fit the Cox/AFT model; please provide ",
             " this dataset in the 'data_Surv' argument of jm().")
    }

    # if the longitudinal outcomes are not in dataS, we set a random value for
    # them. This is needed for the calculation of the matrix of interaction terms
    # between the longitudinal outcomes and other variables.
    for (i in seq_along(respVars)) {
        if (is.null(dataS[[respVars[i]]])) dataS[[respVars[i]]] <- rnorm(nrow(dataS))
    }
    # terms for survival model
    terms_Surv <- Surv_object$terms
    terms_Surv_noResp <- delete.response(terms_Surv)
    mf_surv_dataS <- model.frame.default(terms_Surv, data = dataS)

    # survival times
    Surv_Response <- model.response(mf_surv_dataS)
    type_censoring <- attr(Surv_Response, "type")
    if (is.null(dataS[[idVar]])) {
        if (is.null(id_var)) {
            stop("cannot extract the subject id variable from the dataset used to fit the ",
                 "survival model. Please specify the 'id_var' argument.\n")
        } else {
            idT <- dataS[[id_var]]
        }
    } else {
        idT <- dataS[[idVar]]
    }
    if (!is.null(NAs_surv <- attr(mf_surv_dataS, "na.action"))) {
        idT <- idT[-NAs_surv]
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
    if (!all(order(unique(idT)) == order(unique(dataL[[idVar]])))) {
        warning("It seems that the ordering of the subjects in dataset used to fit the ",
                "mixed models and the dataset used for the survival model is not the same. ",
                "We set internally the datasets in the same order, but it would be best ",
                "that you do it beforehand on your own.")
        dataS <- dataS[order(idT), ]
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
        Time_left <- Time_start <- trunc_Time <- rep(0.0, nT)
        delta <-  unname(Surv_Response[, "status"])
    } else if (type_censoring == "counting") {
        Time_start <- unname(Surv_Response[, "start"])
        Time_stop <- unname(Surv_Response[, "stop"])
        delta <-  unname(Surv_Response[, "status"])
        Time_right <- tapply(Time_stop, idT, tail, n = 1) # time of event
        trunc_Time <- tapply(Time_start, idT, head, n = 1) # possible left truncation time
        Time_left <- rep(0.0, nT)
        delta <- tapply(event, idT, tail, n = 1) # event indicator at Time_right
    } else if (type_censoring == "interval") {
        Time1 <-  unname(Surv_Response[, "time1"])
        Time2 <-  unname(Surv_Response[, "time2"])
        trunc_Time <- Time_start <- rep(0.0, nT)
        delta <- unname(Surv_Response[, "status"])
        Time_right <- Time1
        Time_right[delta == 3] <- Time2[delta == 3]
        Time_right[delta == 2] <- 0.0
        Time_left <- Time1
        Time_left[delta <= 1] <- 0.0
    }
    which_event <- which(delta == 1)
    which_right <- which(delta == 0)
    which_left <- which(delta == 2)
    which_interval <- which(delta == 3)
    # 'Time_integration' is the upper limit of the integral in likelihood
    # of the survival model. For subjects with event (delta = 1), for subjects with
    # right censoring and for subjects with interval censoring we need to integrate
    # up to 'Time_right'. For subjects with left censoring we need to integrate up to
    # 'Time_left'; hence we set for them 'Time_integration = Time_left'. For subjects
    # with interval censoring we need also 'Time_integration2' which is equal to
    # 'Time_left'
    Time_integration <- Time_right
    Time_integration[which_left] <- Time_left[which_left]
    Time_integration2 <- rep(0.0, nT)
    if (length(which_interval)) {
        Time_integration2[which_interval] <- Time_left[which_interval]
    }

    # create Gauss Kronrod points and weights
    GK <- gaussKronrod(con$GK_k)
    sk <- GK$sk
    P <- c(Time_integration - trunc_Time) / 2
    st <- outer(P, sk) + (c(Time_integration + trunc_Time) / 2)
    log_Pwk <- rep(log(P), each = length(sk)) + rep_len(log(GK$wk), length.out = length(st))
    if (length(which_interval)) {
        P2 <- c(Time_integration2 - trunc_Time) / 2
        st2 <- outer(P2, sk) + (c(Time_integration2 + trunc_Time) / 2)
        log_Pwk2 <- rep(log(P2), each = length(sk)) +
            rep_len(log(GK$wk), length.out = length(st2))
    } else {
        P2 <- st2 <- log_Pwk2 <- NULL
    }

    # knots for the log baseline hazard function
    if (is.null(con$knots)) {
        con$knots <- knots(0, floor(max(Time_integration)) + 1,
                           con$base_hazard_segments, con$Bsplines_degree)
    }

    # Extract functional forms per longitudinal outcome
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
    functional_forms <-  functional_forms[order(match(names(functional_forms),
                                                      respVars_form))]
    ###################################################################
    # List of lists
    # One list component per association structure per ouctome
    # List components vectors of integers corresponding to the term
    # each association structure corresponds to
    functional_forms_per_outcome <- lapply(functional_forms,
                                           extract_functional_forms_per_outcome)
    functional_forms_per_outcome <- lapply(functional_forms_per_outcome,
                                           function (x) x[sapply(x, length) > 0])
    collapsed_functional_forms <- lapply(functional_forms_per_outcome, names)
    # NEW INPUT OBJECTS FOR create_Wlong.cpp
    ns_functional_forms_per_outcome <- sapply(functional_forms_per_outcome, length)
    maxs_functional_forms_per_outcome <- sapply(lapply(functional_forms_per_outcome, unlist), max)

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

    W0_H <- splineDesign(con$knots, c(t(st)), ord = con$Bsplines_degree + 1)
    dataS_H <- SurvData_HazardModel(st, dataS, Time_start, idT)
    mf <- model.frame.default(terms_Surv_noResp, data = dataS_H)
    W_H <- model.matrix.default(terms_Surv_noResp, mf)[, -1, drop = FALSE]
    if (!ncol(W_H)) {
        W_H <- matrix(0.0, nrow = nrow(W_H), ncol = 1L)
    }
    X_H <- desing_matrices_functional_forms(st, terms_FE_noResp,
                                            dataL, time_var, idVar,
                                            collapsed_functional_forms)
    Z_H <- desing_matrices_functional_forms(st, terms_RE,
                                            dataL, time_var, idVar,
                                            collapsed_functional_forms)
    U_H <- lapply(seq_along(Mixed_objects), function (i) {
        tt <- terms(functional_forms[[i]])
        model.matrix(tt, model.frame(tt, data = dataS_H))[, -1, drop = FALSE]
    })
    if (length(which_event)) {
        W0_h <- splineDesign(con$knots, Time_right, ord = con$Bsplines_degree + 1)
        dataS_h <- SurvData_HazardModel(Time_right, dataS, Time_start, idT)
        mf <- model.frame.default(terms_Surv_noResp, data = dataS_h)
        W_h <- model.matrix.default(terms_Surv_noResp, mf)[, -1, drop = FALSE]
        if (!ncol(W_h)) {
            W_h <- matrix(0.0, nrow = nrow(W_h), ncol = 1L)
        }
        X_h <- desing_matrices_functional_forms(Time_right, terms_FE_noResp,
                                                dataL, time_var, idVar,
                                                collapsed_functional_forms)
        Z_h <- desing_matrices_functional_forms(Time_right, terms_RE,
                                                dataL, time_var, idVar,
                                                collapsed_functional_forms)
        U_h <- lapply(seq_along(Mixed_objects), function (i) {
            tt <- terms(functional_forms[[i]])
            model.matrix(tt, model.frame(tt, data = dataS_h))[, -1, drop = FALSE]
        })
    } else {
        W0_h <- W_h <- matrix(0.0)
        X_h <- Z_h <- U_h <- rep(list(matrix(0.0)), length(respVars))
    }
    if (length(which_interval)) {
        W0_H2 <- splineDesign(con$knots, c(t(st2)), ord = con$Bsplines_degree + 1)
        dataS_H2 <- SurvData_HazardModel(st2, dataS, Time_start, idT)
        mf2 <- model.frame.default(terms_Surv_noResp, data = dataS_H2)
        W_H2 <- model.matrix.default(terms_Surv_noResp, mf2)[, -1, drop = FALSE]
        if (!ncol(W_H2)) {
            W_H2 <- matrix(0.0, nrow = nrow(W_H2), ncol = 1L)
        }
        X_H2 <- desing_matrices_functional_forms(st, terms_FE_noResp,
                                                 dataL, time_var, idVar,
                                                 collapsed_functional_forms)
        Z_H2 <- desing_matrices_functional_forms(st, terms_RE,
                                                 dataL, time_var, idVar,
                                                 collapsed_functional_forms)
        U_H2 <- lapply(seq_along(Mixed_objects), function (i) {
            tt <- terms(functional_forms[[i]])
            model.matrix(tt, model.frame(tt, data = dataS_H2))[, -1, drop = FALSE]
        })
    } else {
        W0_H2 <- W_H2 <- matrix(0.0)
        X_H2 <- Z_H2 <- U_H2 <- rep(list(matrix(0.0)), length(respVars))
    }
    ######################################################################################
    ######################################################################################
    Data <- list(idL = idL, idL_lp = idL_lp, unq_idL = unq_idL,
                 y = y, X = X, Z = Z, Xhc = Xhc,
                 columns_HC = columns_HC, columns_nHC = columns_nHC,
                 #####
                 idT = idT,
                 Time_right = Time_right, Time_left = Time_left, Time_start = Time_start,
                 delta = delta, which_event = which_event, which_right = which_right,
                 which_left = which_left, which_interval = which_interval,
                 W0_H = W0_H, W_H = W_H, X_H = X_H, Z_H = Z_H, U_H = U_H,
                 W0_h = W0_h, W_h = W_h, X_h = X_h, Z_h = Z_h, U_h = U_h,
                 W0_H2 = W0_H2, W_H2 = W_H2, X_H2 = X_H2, Z_H2 = Z_H2, U_H2 = U_H2,
                 log_Pwk = log_Pwk, log_Pwk2 = log_Pwk2)
    # drop names and other attributes from model matrices
    Data[] <- lapply(Data, drop_names)
    ######################################################################################
    ######################################################################################
    # objects to export
    data <- list(dataL = dataL, dataS = dataS)
    model_info <- list(
        terms = list(terms_FE = terms_FE, terms_FE_noResp = terms_FE_noResp,
                     terms_RE = terms_RE, terms_Surv,
                     terms_Surv_noResp = terms_Surv_noResp),
        var_names = list(respVars = respVars, respVars_form = respVars_form,
                         idVar = idVar, time_var = time_var),
        families = families,
        ids = list(idL = idL, idL_lp = idL_lp, unq_idL = unq_idL, idT = idT),
        n = nY,
        HC = list(columns_HC = columns_HC, columns_nHC = columns_nHC),
        type_censoring = type_censoring,
        fun_forms = list(functional_forms = functional_forms,
                         functional_forms_per_outcome = functional_forms_per_outcome,
                         collapsed_functional_forms = collapsed_functional_forms)
    )
    ######################################################################################
    ######################################################################################
    # initial values
    betas <- lapply(Mixed_objects, fixef)
    log_sigmas <- lapply(Mixed_objects, extract_log_sigmas)
    D_lis <- lapply(Mixed_objects, extract_D)
    D <- bdiag(D_lis)
    b <- mapply(extract_b, Mixed_objects, unq_idL, MoreArgs = list(n = nY),
                SIMPLIFY = FALSE)
    init_surv <- init_vals_surv(Data, model_info, data, betas, b, con)
    bs_gammas <- init_surv$bs_gammas
    gammas <- init_surv$gammas
    alphas <- init_surv$alphas
    # we are going to have multiple chains; hence, we need to randomly permute these
    # initial values
    initial_values <- list(betas = betas, log_sigmas = log_sigmas, D = D,
                           b = b, bs_gammas = bs_gammas, gammas = gammas,
                           alphas = alphas)
    ######################################################################################
    ######################################################################################
    # variance covariance matrices for proposal distributions in
    # the Metropolis-Hastings algorithm
    #  - betas the fixed effects that in the hierarchical centering part
    #  - tilde_betas the fixed effects that are not in the hierarchical centering part
    vcov_prop_betas <- mapply(get_vcov_FE, Mixed_objects, columns_nHC,
                              MoreArgs = list(which = "betas"), SIMPLIFY = FALSE)
    vcov_prop_tilde_betas <- mapply(get_vcov_FE, Mixed_objects, columns_nHC,
                                    MoreArgs = list(which = "tilde_betas"),
                                    SIMPLIFY = FALSE)
    r <- mapply(extract_vcov_prop_RE, Mixed_objects, Z, idL, SIMPLIFY = FALSE)
    vcov_prop_RE <- array(0.0, c(dim(D), nY))
    for (i in seq_len(nY)) {
        rr <- lapply(r, function (m, i) m[[as.character(i)]], i = i)
        if (any(ind <- sapply(rr, is.null))) rr[ind] <- D_lis[ind]
        vcov_prop_RE[, , i] <- .bdiag(rr)
    }
    vcov_prop_bs_gammas <- init_surv$vcov_prop_bs_gammas
    vcov_prop_gammas <- init_surv$vcov_prop_gammas
    vcov_prop_alphas <- init_surv$vcov_prop_alphas
    vcov_prop <- list(vcov_prop_betas = vcov_prop_betas,
                      vcov_prop_tilde_betas = vcov_prop_tilde_betas,
                      vcov_prop_RE = vcov_prop_RE,
                      vcov_prop_bs_gammas = vcov_prop_bs_gammas,
                      vcov_prop_gammas = vcov_prop_gammas,
                      vcov_prop_alphas = vcov_prop_alphas)
    list(betas = betas, D = D, log_sigmas = log_sigmas,
         bs_gammas = bs_gammas, gammas = gammas, alphas = alphas,
         vcov_prop = vcov_prop)
}

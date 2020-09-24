jm <- function (Surv_object, Mixed_objects, time_var,
                functional_forms = NULL, data_Surv = NULL, id_var = NULL,
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
    # - cores: the number of cores to use for running the chains in parallel
    # - MALA: if TRUE, the MALA algorithm is used when update the elements of
    #         of the Cholesky factor of the D matrix
    con <- list(GK_k = 15L, Bsplines_degree = 2, base_hazard_segments = 10,
                diff = 2L, n_chains = 3L, n_burnin = 500L, n_iter = 3500L,
                seed = 123L,  cores = max(parallel::detectCores() - 1, 1),
                MALA = FALSE)
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
    idL_ind <- lapply(idL, function (x) seq_along(x))
    idL_ind <- mapply2(function (x, y) split(x, y), idL_ind, idL)
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
    mf_FE_dataL <- mapply(fix_NAs_fixed, mf_FE_dataL, NAs_FE_dataL, NAs_RE_dataL,
                          SIMPLIFY = FALSE)
    mf_RE_dataL <- mapply(fix_NAs_random, mf_RE_dataL, NAs_RE_dataL, NAs_FE_dataL,
                          SIMPLIFY = FALSE)

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
    X <- mapply(model.matrix.default, terms_FE, mf_FE_dataL, SIMPLIFY = FALSE)
    Z <- mapply(model.matrix.default, terms_RE, mf_RE_dataL, SIMPLIFY = FALSE)
    if (length(Z) == 1 && ncol(Z[[1]]) == 1) {
        stop("jm() does not currently work when you have a single ",
             "longitudinal outcome and only random intercepts.")
    }
    ind_RE <- lapply(Z, FUN = function(x) seq_len(ncol(x)))
    if (length(ind_RE) > 1) {
        for (i in 2:length(ind_RE)) {
            ind_RE[[i]] <- max(ind_RE[[i - 1]])
        }
    }
    componentsHC <- mapply2(create_HC_X2, X, Z, idL)
    Xbase <- lapply(componentsHC, "[[", "Xbase")
    baseline <- lapply(componentsHC, "[[", "baseline")
    x_in_z <- lapply(componentsHC, "[[", "x_in_z")
    x_notin_z <- lapply(componentsHC, "[[", "x_notin_z")
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
                 "survival model. Please include this variable in the dataset ",
                 "and/or specify the 'id_var' argument.\n")
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
        warning("It seems that the ordering of the subjects in the dataset used to fit the ",
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
        delta <- tapply(delta, idT, tail, n = 1) # event indicator at Time_right
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
    Time_integration2 <- rep(0.0, nT)
    if (length(which_interval)) {
        Time_integration2[which_interval] <- Time_right[which_interval]
    }

    # create Gauss Kronrod points and weights
    GK <- gaussKronrod(con$GK_k)
    sk <- GK$sk
    P <- c(Time_integration - trunc_Time) / 2
    st <- outer(P, sk) + (c(Time_integration + trunc_Time) / 2)
    log_Pwk <- rep(log(P), each = length(sk)) + rep_len(log(GK$wk), length.out = length(st))
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
        qs <- quantile(c(Time_right, Time_left), probs = c(0.1, 0.9))
        con$knots <- knots(qs[1L], qs[2L], con$base_hazard_segments,
                           con$Bsplines_degree)
    }
    .knots_baseline_hazard <- con$knots
    env <- new.env(parent = .GlobalEnv)
    assign(".knots_baseline_hazard", con$knots, envir = env)

    # create long version of idT
    idTs <- rep(idT, each = con$GK_k)

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
    functional_forms <- functional_forms[order(match(names(functional_forms),
                                                     respVars_form))]
    ###################################################################
    # List of lists
    # One list component per association structure per outcome
    # List components vectors of integers corresponding to the term
    # each association structure corresponds to
    FunForms_per_outcome <- lapply(functional_forms, extract_functional_forms,
                                   data = dataS)
    FunForms_per_outcome <- lapply(FunForms_per_outcome,
                                   function (x) x[sapply(x, length) > 0])
    collapsed_functional_forms <- lapply(FunForms_per_outcome, names)

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
    W0_H <- splineDesign(con$knots, c(t(st)), ord = con$Bsplines_degree + 1,
                         outer.ok = TRUE)
    dataS_H <- SurvData_HazardModel(st, dataS, Time_start, idT, time_var)
    mf <- model.frame.default(terms_Surv_noResp, data = dataS_H)
    W_H <- model.matrix.default(terms_Surv_noResp, mf)[, -1, drop = FALSE]
    any_gammas <- as.logical(ncol(W_H))
    if (!any_gammas) {
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
        W0_h <- splineDesign(con$knots, Time_right,
                             ord = con$Bsplines_degree + 1, outer.ok = TRUE)
        dataS_h <- SurvData_HazardModel(Time_right, dataS, Time_start, idT,
                                        time_var)
        mf <- model.frame.default(terms_Surv_noResp, data = dataS_h)
        W_h <- model.matrix.default(terms_Surv_noResp, mf)[, -1, drop = FALSE]
        if (!any_gammas) {
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
        W0_H2 <- splineDesign(con$knots, c(t(st2)),
                              ord = con$Bsplines_degree + 1, outer.ok = TRUE)
        dataS_H2 <- SurvData_HazardModel(st2, dataS, Time_start, idT, time_var)
        mf2 <- model.frame.default(terms_Surv_noResp, data = dataS_H2)
        W_H2 <- model.matrix.default(terms_Surv_noResp, mf2)[, -1, drop = FALSE]
        if (!any_gammas) {
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
    ############################################################################
    ############################################################################
    Data <- list(n = nY, idL = idL, idL_ind = idL_ind, idL_lp = idL_lp, unq_idL = unq_idL,
                 y = y, X = X, Z = Z, Xbase = Xbase,
                 baseline = baseline, x_in_z = x_in_z, x_notin_z = x_notin_z,
                 #####
                 idT = idT, any_gammas = any_gammas,
                 Time_right = Time_right, Time_left = Time_left, Time_start = Time_start,
                 delta = delta, which_event = which_event, which_right = which_right,
                 which_left = which_left, which_interval = which_interval,
                 W0_H = W0_H, W_H = W_H, X_H = X_H, Z_H = Z_H, U_H = U_H,
                 W0_h = W0_h, W_h = W_h, X_h = X_h, Z_h = Z_h, U_h = U_h,
                 W0_H2 = W0_H2, W_H2 = W_H2, X_H2 = X_H2, Z_H2 = Z_H2, U_H2 = U_H2,
                 log_Pwk = log_Pwk, log_Pwk2 = log_Pwk2,
                 ind_RE = ind_RE, extra_parms = rep(0.0, length(y)))
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
                         idVar = idVar, time_var = time_var),
        families = families,
        type_censoring = type_censoring,
        functional_forms = functional_forms,
        FunForms_per_outcome = FunForms_per_outcome,
        collapsed_functional_forms = collapsed_functional_forms,
        FunForms_cpp = lapply(FunForms_per_outcome, unlist),
        FunForms_ind = FunForms_ind(FunForms_per_outcome)
    )
    ############################################################################
    ############################################################################
    # initial values
    betas <- lapply(Mixed_objects, fixef)
    log_sigmas <- sapply(Mixed_objects, extract_log_sigmas)
    D_lis <- lapply(Mixed_objects, extract_D)
    D <- bdiag(D_lis)
    b <- mapply2(extract_b, Mixed_objects, unq_idL, MoreArgs = list(n = nY))
    init_surv <- init_vals_surv(Data, model_info, data, betas, b, con)
    bs_gammas <- init_surv$bs_gammas
    gammas <- init_surv$gammas
    alphas <- init_surv$alphas
    initial_values <- list(betas = betas, log_sigmas = log_sigmas, D = D,
                           b = b, bs_gammas = bs_gammas, gammas = gammas,
                           alphas = alphas, tau_bs_gammas = 20)
    ############################################################################
    ############################################################################
    # variance covariance matrices for proposal distributions in
    # the Metropolis-Hastings algorithm
    #  - betas the fixed effects that in the hierarchical centering part
    #  - tilde_betas the fixed effects that are not in the hierarchical
    #    centering part
    vcov_prop_betas <- lapply(Mixed_objects, vcov2)
    r <- mapply2(extract_vcov_prop_RE, Mixed_objects, Z, idL)
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
                      vcov_prop_RE = vcov_prop_RE,
                      vcov_prop_bs_gammas = vcov_prop_bs_gammas,
                      vcov_prop_gammas = vcov_prop_gammas,
                      vcov_prop_alphas = vcov_prop_alphas)
    ############################################################################
    ############################################################################
    # Priors
    DD <- diag(ncol(W0_H))
    Tau_bs_gammas <- crossprod(diff(DD, differences = con$diff))
    prs <- list(mean_betas = lapply(betas, "*", 0.0),
                Tau_betas = lapply(betas, function (b) 0.01 * diag(length(b))),
                mean_gammas = gammas * 0.0,
                Tau_gammas = 0.01 * diag(length(gammas)),
                mean_bs_gammas = bs_gammas * 0.0,
                Tau_bs_gammas = Tau_bs_gammas,
                mean_alphas = lapply(alphas, "*", 0.0),
                Tau_alphas = lapply(alphas, function (a) 0.01 * diag(length(a))),
                A_tau_bs_gammas = 1, B_tau_bs_gammas = 0.1,
                rank_Tau_bs_gammas = qr(Tau_bs_gammas)$rank,
                prior_D_sds_df = 3.0, prior_D_sds_sigma = 10.0,
                prior_D_L_etaLKJ = 3.0, prior_sigmas_df = 3.0,
                prior_sigmas_sigma = 20.0)
    if (is.null(priors) || !is.list(priors)) {
        priors <- prs
    } else {
        ind <- names(prs) %in% names(priors)
        prs[ind] <- priors[ind]
        priors <- prs
    }
    ############################################################################
    ############################################################################
    # Fit the model
    out <- jm_fit(Data, model_info, initial_values, priors, con, vcov_prop)
    S <- lapply(out$mcmc, summary)
    statistics <- list(
        Mean = lapply(S, get_statistic, "Mean"),
        Median = lapply(S, get_statistic, "Median"),
        SD = lapply(S, get_statistic, "SD"),
        SE = lapply(S, get_statistic, "Time-series SE"),
        CI_low = lapply(S, get_statistic, "2.5CI"),
        CI_upp = lapply(S, get_statistic, "97.5CI"),
        P = lapply(out$mcmc, function (x) apply(do.call("rbind", x), 2L, Ptail)),
        Effective_Size = lapply(out$mcmc, function (x)
            apply(do.call("rbind", x), 2L, effective_size))
    )
    if (!is.null(out$mcmc[["b"]])) {
        znams <- unlist(lapply(Data$Z, colnames), use.names = FALSE)
        l <- sapply(Data$unq_idL, length)
        dnames_b <- list(unlist(Data$unq_idL[which.max(l)]), znams)
        fix_b <- function (stats) {
            x <- stats$b
            dim(x) <- sapply(dnames_b, length)
            dimnames(x) <- dnames_b
            stats$b <- x
            stats
        }
        statistics[] <- lapply(statistics, fix_b)
    }
    if (con$n_chains > 1) {
        statistics <- c(statistics,
                        Rhat = list(lapply(out$mcm, function (theta)
                            coda::gelman.diag(theta)$psrf)))
    }
    out <- c(out, list(statistics = statistics,
                       model_data = Data, model_info = model_info,
                       initial_values = initial_values,
                       control = con, priors = priors, call = call,
                       vcov_prop = vcov_prop))
    class(out) <- "jm"
    out
}

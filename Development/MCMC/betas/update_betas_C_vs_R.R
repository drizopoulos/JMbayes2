#setwd("V:/Users/051528(PM_Afonso)/GitHub/JMbayes2")
library("survival")
library("nlme")
library("GLMMadaptive")
library("coda")
library("lattice")
library("splines")
source("./R/jm.R")
source("./R/jm_fit.R")
source("./R/help_functions.R")
source("./R/basic_methods.R")
source("./R/create_Wlong_mats.R")
source("./Development/jm/PBC_data.R")

pbc2$prothrombin[pbc2$id == levels(pbc2$id)[1L]] <- NA
pbc2$serBilir[pbc2$id == levels(pbc2$id)[1L]] <- NA
pbc2$ascites[pbc2$id == levels(pbc2$id)[1L]] <- NA

pbc2$prothrombin[pbc2$id == levels(pbc2$id)[2L]] <- NA

fm1 <- lme(log(serBilir) ~ year * (drug + sex) + I(year^2) + age + serChol,
           data = pbc2, random = ~ year + I(year^2)| id, na.action = na.exclude)

fm2 <- lme(prothrombin ~ ns(year, 2) + sex, data = pbc2,
           random = ~ year + I(year^2)| id,
           na.action = na.exclude, control = lmeControl(opt = "optim"))

fm3 <- mixed_model(ascites ~ year + sex, data = pbc2, random = ~ year | id,
                   family = binomial())

Mixed <- list(fm1, fm2, fm3)
Cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

################################################################################

Surv_object = Cox
Mixed_objects = Mixed
time_var = 'year'
functional_forms = NULL
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL

# jm()
{
  con <- list(GK_k = 15L, Bsplines_degree = 2, base_hazard_segments = 10,
              diff = 2L, n_chains = 3L, n_burnin = 500L, n_iter = 3500L,
              n_thin = 1L, seed = 123L, MALA = FALSE,
              save_random_effects = FALSE,
              cores = max(parallel::detectCores() - 1, 1))
  #control <- c(control, list(...))
  namC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% namC]) > 0)
    warning("unknown names in control: ", paste(noNms, collapse = ", "))
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
  idL_ind <- lapply(idL, function (x) seq_along(x))
  idL_ind <- mapply2(function (x, y) split(x, y), idL_ind, idL)
  nY <- length(unique(idL))
  # order data by idL and time_var
  if (is.null(dataL[[time_var]])) {
    stop("the variable specified in agument 'time_var' cannot be found ",
         "in the database of the longitudinal models.")
  }
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
  nres <- sapply(Z, ncol)
  ind_RE <- split(seq_len(sum(nres)), rep(seq_along(Z), nres))
  componentsHC <- mapply2(create_HC_X2, X, Z, idL)
  Xbase <- lapply(componentsHC, "[[", "Xbase")
  #Xbase[] <- mapply2(function (m, nams) {rownames(m) <- nams; m}, Xbase, unq_idL)
  baseline <- lapply(componentsHC, "[[", "baseline")
  x_in_z <- lapply(componentsHC, "[[", "x_in_z")
  x_notin_z <- lapply(componentsHC, "[[", "x_notin_z")
  nfes <- sapply(X, ncol)
  # 'ind_FE' is used in vec2field() to re-create the field of betas
  # from betas_vec
  ind_FE <- split(seq_len(sum(nfes)), rep(seq_along(X), nfes))
  x_in_z_base <- mapply2(function (x, y) sort(c(x, y)), x_in_z, baseline)
  # 'ind_FE_HC' denotes which elements of betas_vec are in the HC formulation
  # this will be use to save the results in the corresponding columns
  ind_FE_HC <- unlist(mapply2(function (x, ind) x[ind], ind_FE, x_in_z_base),
                      use.names = FALSE)
  # 'ind_FE_HC' denotes which elements of betas_vec are not in the
  # HC formulation. It is a list, to be used in conjuction with
  # has_tilde_betas. That is, if has_tilde_betas = TRUE, we need to save
  # from the Metropolis-Hastings step the betas in the columns 'ind_FE_nHC'
  ind_FE_nHC <- mapply2(function (x, ind) x[-ind], ind_FE, x_in_z_base)
  has_tilde_betas <- as.integer(sapply(ind_FE_nHC, length) > 0)
  ind_FE_nHC[] <- lapply(ind_FE_nHC, function (x) if (length(x)) x else 0L)
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
  # if the time_var is not in dataS set it to a random number
  if (is.null(dataS[[time_var]])) dataS[[time_var]] <- rnorm(nrow(dataS))
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
    dataS <- dataS[-NAs_surv, ]
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
    Time_left <- Time_start <- trunc_Time <- rep(0.0, nrow(dataS))
    delta <-  unname(Surv_Response[, "status"])
  } else if (type_censoring == "counting") {
    Time_start <- unname(Surv_Response[, "start"])
    Time_stop <- unname(Surv_Response[, "stop"])
    delta <-  unname(Surv_Response[, "status"])
    Time_right <- tapply(Time_stop, idT, tail, n = 1) # time of event
    trunc_Time <- tapply(Time_start, idT, head, n = 1) # possible left truncation time
    Time_left <- rep(0.0, nrow(dataS))
    delta <- tapply(delta, idT, tail, n = 1) # event indicator at Time_right
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
  names(Time_right) <- names(Time_left) <- names(Time_start) <- idT
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
  n_strata <- length(unique(strata))
  
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
  .knots_base_hazard <- con$knots
  env <- new.env(parent = .GlobalEnv)
  assign(".knots_base_hazard", con$knots, envir = env)
  
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
  set_env <- function (form, env) {environment(form) <- env; form}
  functional_forms[] <- lapply(functional_forms, set_env, env = env)
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
  strata_H <- rep(strata, each = con$GK_k)
  W0_H <- create_W0(c(t(st)), con$knots, con$Bsplines_degree + 1, strata_H)
  dataS_H <- SurvData_HazardModel(st, dataS, Time_start,
                                  paste0(idT, "_", strata), time_var)
  mf <- model.frame.default(terms_Surv_noResp, data = dataS_H)
  W_H <- construct_Wmat(terms_Surv_noResp, mf)
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
  U_H <- lapply(functional_forms, construct_Umat, dataS = dataS_H)
  if (length(which_event)) {
    W0_h <- create_W0(Time_right, con$knots, con$Bsplines_degree + 1, strata)
    dataS_h <- SurvData_HazardModel(Time_right, dataS, Time_start,
                                    paste0(idT, "_", strata), time_var)
    mf <- model.frame.default(terms_Surv_noResp, data = dataS_h)
    W_h <- construct_Wmat(terms_Surv_noResp, mf)
    if (!any_gammas) {
      W_h <- matrix(0.0, nrow = nrow(W_h), ncol = 1L)
    }
    X_h <- desing_matrices_functional_forms(Time_right, terms_FE_noResp,
                                            dataL, time_var, idVar,
                                            collapsed_functional_forms)
    Z_h <- desing_matrices_functional_forms(Time_right, terms_RE,
                                            dataL, time_var, idVar,
                                            collapsed_functional_forms)
    U_h <- lapply(functional_forms, construct_Umat, dataS = dataS_h)
  } else {
    W0_h <- W_h <- matrix(0.0)
    X_h <- Z_h <- U_h <- rep(list(matrix(0.0)), length(respVars))
  }
  if (length(which_interval)) {
    W0_H2 <- create_W0(c(t(st2)), con$knots, con$Bsplines_degree + 1,
                       strata_H)
    dataS_H2 <- SurvData_HazardModel(st2, dataS, Time_start,
                                     paste0(idT, "_", strata), time_var)
    mf2 <- model.frame.default(terms_Surv_noResp, data = dataS_H2)
    W_h <- construct_Wmat(terms_Surv_noResp, mf2)
    if (!any_gammas) {
      W_H2 <- matrix(0.0, nrow = nrow(W_H2), ncol = 1L)
    }
    X_H2 <- desing_matrices_functional_forms(st, terms_FE_noResp,
                                             dataL, time_var, idVar,
                                             collapsed_functional_forms)
    Z_H2 <- desing_matrices_functional_forms(st, terms_RE,
                                             dataL, time_var, idVar,
                                             collapsed_functional_forms)
    U_H <- lapply(functional_forms, construct_Umat, dataS = dataS_H2)
  } else {
    W0_H2 <- W_H2 <- matrix(0.0)
    X_H2 <- Z_H2 <- U_H2 <- rep(list(matrix(0.0)), length(respVars))
  }
  nfes_HC <- sapply(x_in_z_base, length)
  out_in <- sapply(idL, "%in%", x = seq_len(nT))
  all_pat <- apply(out_in, 1L, paste0, collapse = "/")
  id_patt <- match(all_pat, unique(all_pat))
  find_patt <- function (patt, n) which(rep(patt, times = n))
  ind_RE_patt <- apply(unique(out_in), 1L, find_patt, n = nres)
  ind_FE_patt <- apply(unique(out_in), 1L, find_patt, n = nfes_HC)
  # X_dot <- create_X_dot2(nT, nres, ind_FE_HC, x_in_z, x_in_z_base, unq_idL,
  #                        Xbase)
  X_dot <- create_X_dot(Xbase, nT, unq_idL, nres, nfes_HC, baseline, x_in_z_base, x_in_z)
  ############################################################################
  ############################################################################
  Data <- list(n = nY, idL = idL, idL_ind = idL_ind, idL_lp = idL_lp, unq_idL = unq_idL,
               y = y, X = X, Z = Z, Xbase = Xbase, baseline = baseline, X_dot = X_dot,
               x_in_z = x_in_z, x_notin_z = x_notin_z,
               has_tilde_betas = has_tilde_betas, ind_FE = ind_FE,
               ind_FE_HC = ind_FE_HC, ind_FE_nHC = ind_FE_nHC,
               id_patt = id_patt, ind_RE_patt = ind_RE_patt,
               ind_FE_patt = ind_FE_patt,
               #####
               idT = idT, any_gammas = any_gammas, strata = strata,
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
                         alphas = alphas, tau_bs_gammas = rep(20, n_strata))
  ############################################################################
  ############################################################################
  # variance covariance matrices for proposal distributions in
  # the Metropolis-Hastings algorithm
  #  - betas the fixed effects that in the hierarchical centering part
  #  - tilde_betas the fixed effects that are not in the hierarchical
  #    centering part
  vcov_prop_betas <- lapply(Mixed_objects, vcov2)
  vcov_prop_betas <- mapply2(get_betas_nHC, vcov_prop_betas, x_notin_z)
  r <- mapply2(extract_vcov_prop_RE, Mixed_objects, Z, idL)
  vcov_prop_RE <- array(0.0, c(dim(D), nY))
  for (i in seq_len(nY)) {
    rr <- lapply(r, function (m, i) m[[as.character(i)]], i = i)
    if (any(ind <- sapply(rr, is.null))) rr[ind] <- lapply(D_lis[ind], "*", 0.1)
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
  # we define as thetas = c(gammas, unlist(alphas))
  # for some of these coefficients we place penalties/shrinkage priors
  # 'Tau_thetas' is list of the prior precision matrices for thetas
  # 'ind_thetas' is a list of position indicators specifying which of the
  # thetas are involved in the calculation of the corresponding tau parameter,
  # where tau is the penalty parameter. Hence, the number of taus we will need
  # to estimate is the length of ind_thetas.
  # For each of these taus we place a Gamma prior with parameters A_tau and
  # B_tau. Thus, A_tau and B_tau are vectors with length equal to ind_thetas
  #
  thetas <- if (any_gammas) {
    c(gammas, unlist(alphas))
  } else {
    c(unlist(alphas))
  }
  A_tau <- rep(1, n_strata)
  B_tau <- rep(0.1, n_strata)
  ind_thetas <- list(alphas = grep("tv(", unlist(lapply(U_H, colnames)),
                                   fixed = TRUE))
  Tau_bs_gammas <- crossprod(diff(diag(ncol(W0_H) / n_strata),
                                  differences = con$diff))
  Tau_bs_gammas <- rep(list(Tau_bs_gammas), n_strata)
  mean_betas <- lapply(betas, "*", 0.0)
  Tau_betas <- lapply(betas, function (b) 0.01 * diag(length(b)))
  mean_betas_HC <- unlist(mean_betas, use.names = FALSE)[ind_FE_HC]
  Tau_betas_HC <- bdiag(Tau_betas)[ind_FE_HC, ind_FE_HC, drop = FALSE]
  mean_betas_nHC <- mapply2(get_betas_nHC, mean_betas, x_notin_z)
  Tau_betas_nHC <- mapply2(get_betas_nHC, Tau_betas, x_notin_z)
  prs <- list(mean_betas_HC = mean_betas_HC, Tau_betas_HC = Tau_betas_HC,
              mean_betas_nHC = mean_betas_nHC, Tau_betas_nHC = Tau_betas_nHC,
              mean_bs_gammas = lapply(Tau_bs_gammas, function (x) x[, 1] * 0),
              Tau_bs_gammas = Tau_bs_gammas,
              A_tau_bs_gammas = rep(1, n_strata), B_tau_bs_gammas = rep(0.1, n_strata),
              rank_Tau_bs_gammas =
                sapply(lapply(Tau_bs_gammas, qr), "[[", 'rank'),
              mean_gammas = gammas * 0.0, Tau_gammas = 0.1 * diag(length(gammas)),
              penalty_gammas = "none",
              A_lambda_gammas = 0.5, B_lambda_gammas = 1,
              A_tau_gammas = 0.5, B_tau_gammas = 1,
              A_nu_gammas = 0.5, B_nu_gammas = 1,
              A_xi_gammas = 0.5, B_xi_gammas = 1,
              mean_alphas = lapply(alphas, "*", 0.0),
              Tau_alphas = lapply(alphas, function (a) 0.1 * diag(length(a))),
              penalty_alphas = "none",
              A_lambda_alphas = 0.5, B_lambda_alphas = 1,
              A_tau_alphas = 0.5, B_tau_alphas = 1,
              A_nu_alphas = 0.5, B_nu_alphas = 1,
              A_xi_alphas = 0.5, B_xi_alphas = 1,
              prior_D_sds_df = 3.0, prior_D_sds_sigma = 10.0,
              prior_D_L_etaLKJ = 3.0, prior_sigmas_df = 3.0,
              prior_sigmas_sigma = 20.0)
  if (is.null(priors) || !is.list(priors)) {
    priors <- prs
  } else {
    ind <- names(prs) %in% names(priors)
    prs[ind] <- priors
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
}

model_data <- Data
control <- con
control$n_chains <- 1

# jm_fit()
{
  # extract family names
  model_info$family_names <- sapply(model_info$families, "[[", "family")
  # extract link names
  model_info$links <- sapply(model_info$families, "[[", "link")
  # set each y element to a matrix
  model_data$y[] <- lapply(model_data$y, as.matrix)
  # for family = binomial and when y has two columns, set the second column
  # to the number of trials instead the number of failures
  binomial_data <- model_info$family_names == "binomial"
  trials_fun <- function (y) {
    if (NCOL(y) == 2L) y[, 2L] <- y[, 1L] + y[, 2L]
    y
  }
  model_data$y[binomial_data] <-
    lapply(model_data$y[binomial_data], trials_fun)
  idT <- model_data$idT
  id_H <- rep(paste0(idT, "_", model_data$strata), each = control$GK_k)
  id_H <- match(id_H, unique(id_H))
  id_H_ <- rep(idT, each = control$GK_k)
  id_H_ <- match(id_H_, unique(id_H_))
  id_h <- unclass(idT)
  model_data <-
    c(model_data, create_Wlong_mats(model_data, model_info,
                                    initial_values, priors,
                                    control),
      list(id_H = id_H, id_H_ = id_H_, id_h = id_h))
  # cbind the elements of X_H and Z_H, etc.
  model_data$X_H[] <- lapply(model_data$X_H, docall_cbind)
  model_data$X_h[] <- lapply(model_data$X_h, docall_cbind)
  model_data$X_H2[] <- lapply(model_data$X_H2, docall_cbind)
  model_data$Z_H[] <- lapply(model_data$Z_H, docall_cbind)
  model_data$Z_h[] <- lapply(model_data$Z_h, docall_cbind)
  model_data$Z_H2[] <- lapply(model_data$Z_H2, docall_cbind)
  # center the design matrices for the baseline covariates and
  # the longitudinal process
  model_data$W_bar <- rbind(colMeans(model_data$W_H))
  model_data$W_H <- center_fun(model_data$W_H, model_data$W_bar)
  model_data$W_h <- center_fun(model_data$W_h, model_data$W_bar)
  model_data$W_H2 <- center_fun(model_data$W_H2, model_data$W_bar)
  
  model_data$Wlong_bar <- lapply(model_data$Wlong_H, colMeans)
  model_data$Wlong_H <- mapply2(center_fun, model_data$Wlong_H,
                                model_data$Wlong_bar)
  model_data$Wlong_h <- mapply2(center_fun, model_data$Wlong_h,
                                model_data$Wlong_bar)
  model_data$Wlong_H2 <- mapply2(center_fun, model_data$Wlong_H2,
                                 model_data$Wlong_bar)
  model_data$Wlong_bar <- lapply(model_data$Wlong_bar, rbind)
  # unlist priors and initial values for alphas
  initial_values$alphas <- unlist(initial_values$alphas, use.names = FALSE)
  priors$mean_alphas <- unlist(priors$mean_alphas, use.names = FALSE)
  priors$Tau_alphas <- .bdiag(priors$Tau_alphas)
}

################################################################################
# C
control$n_iter <- 1

testC <- update_betas_Gibbs (control, 
                            initial_values,
                            priors,
                            model_data)

# mean_1: 
# 0.2979
# 0.1637
# -0.1079
# -0.2937
# -0.0025
# 0.0035
# 0.0010
# -0.0102
# 0.3511
# 11.2068
# 3.9195
# 
# Sigma_1: 
#   9.9305e-002  1.2307e-004 -2.5370e-003 -3.1478e-002 -3.2242e-006 -1.3526e-003            0            0            0            0            0
# 1.2307e-004  2.2989e-004  3.1255e-009  3.8780e-008 -1.5185e-005  1.6664e-009            0            0            0            0            0
# -2.5370e-003  3.1255e-009  1.0473e-002  4.7715e-004 -8.1879e-011 -6.3549e-005            0            0            0            0            0
# -3.1478e-002  3.8780e-008  4.7715e-004  2.5786e-002 -1.0159e-009  1.7062e-004            0            0            0            0            0
# -3.2242e-006 -1.5185e-005 -8.1879e-011 -1.0159e-009  1.1731e-006 -4.3655e-011            0            0            0            0            0
# -1.3526e-003  1.6664e-009 -6.3549e-005  1.7062e-004 -4.3655e-011  2.4732e-005            0            0            0            0            0
# 0            0            0            0            0            0  1.5057e-002 -1.5057e-002            0            0            0
# 0            0            0            0            0            0 -1.5057e-002  1.7036e-002            0            0            0
# 0            0            0            0            0            0            0            0  1.2629e-001  7.3744e-004 -1.2482e-001
# 0            0            0            0            0            0            0            0  7.3744e-004  3.7442e-004  9.2163e-007
# 0            0            0            0            0            0            0            0 -1.2482e-001  9.2163e-007  1.4118e-001
# 
# D_inv.at(0): [matrix size: 0x0]
# 
# D_inv.at(1): 
#   1.2847e+000 -3.1378e+000 -3.7082e+001            0            0
# -3.1378e+000  1.0667e+002  1.3721e+003            0            0
# -3.7082e+001  1.3721e+003  2.0473e+004            0            0
# 0            0            0  2.2195e-001 -4.3768e-001
# 0            0            0 -4.3768e-001  9.4506e+000
# 
# D_inv.at(2): 
#   1.2847e+000 -3.1378e+000 -3.7082e+001            0            0            0            0            0
# -3.1378e+000  1.0667e+002  1.3721e+003            0            0            0            0            0
# -3.7082e+001  1.3721e+003  2.0473e+004            0            0            0            0            0
# 0            0            0  1.8442e+000  1.7256e-001 -2.2557e+000            0            0
# 0            0            0  1.7256e-001  4.1285e+001 -5.4413e+002            0            0
# 0            0            0 -2.2557e+000 -5.4413e+002  1.0328e+004            0            0
# 0            0            0            0            0            0  2.2195e-001 -4.3768e-001
# 0            0            0            0            0            0 -4.3768e-001  9.4506e+000
# 
# D_inv.at(2): 
#   1.8442e+000  1.7256e-001 -2.2557e+000            0            0
# 1.7256e-001  4.1285e+001 -5.4413e+002            0            0
# -2.2557e+000 -5.4413e+002  1.0328e+004            0            0
# 0            0            0  2.2195e-001 -4.3768e-001
# 0            0            0 -4.3768e-001  9.4506e+000

################################################################################

prior_mean_betas_HC <- mean_betas_HC
prior_Tau_betas_HC <- Tau_betas_HC

b_mat <- do.call(cbind, b)
prior_mean_betas_HC <- mean_betas_HC
prior_Tau_betas_HC <- Tau_betas_HC

testR <- update_betasR(betas, 
                       prior_mean_betas_HC, prior_Tau_betas_HC,
                       b_mat, D, X_dot,
                       ind_FE,
                       ind_FE_HC,
                       id_patt,
                       ind_RE_patt,
                       ind_FE_patt,
                       it = 1,
                       y)

round(testR$mean_1, 4)
testR$Sigma_1
testR$D_inv

# > round(testR$mean_1, 3)
# [,1]
# [1,]   0.2979
# [2,]   0.1637
# [3,]  -0.1079
# [4,]  -0.2937
# [5,]  -0.0025
# [6,]   0.0035
# [7,]   0.0010
# [8,]  -0.0102
# [9,]   0.3511
# [10,]  11.2068
# [11,]  3.919
# 
# > testR$Sigma_1
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]        [,7]        [,8]         [,9]        [,10]         [,11]
# [1,]  9.930536e-02  1.230737e-04 -2.537020e-03 -3.147810e-02 -3.224158e-06 -1.352646e-03  0.00000000  0.00000000  0.000000000 0.000000e+00  0.000000e+00
# [2,]  1.230737e-04  2.298879e-04  3.125507e-09  3.877977e-08 -1.518471e-05  1.666406e-09  0.00000000  0.00000000  0.000000000 0.000000e+00  0.000000e+00
# [3,] -2.537020e-03  3.125507e-09  1.047321e-02  4.771523e-04 -8.187882e-11 -6.354935e-05  0.00000000  0.00000000  0.000000000 0.000000e+00  0.000000e+00
# [4,] -3.147810e-02  3.877977e-08  4.771523e-04  2.578587e-02 -1.015913e-09  1.706200e-04  0.00000000  0.00000000  0.000000000 0.000000e+00  0.000000e+00
# [5,] -3.224158e-06 -1.518471e-05 -8.187882e-11 -1.015913e-09  1.173092e-06 -4.365479e-11  0.00000000  0.00000000  0.000000000 0.000000e+00  0.000000e+00
# [6,] -1.352646e-03  1.666406e-09 -6.354935e-05  1.706200e-04 -4.365479e-11  2.473244e-05  0.00000000  0.00000000  0.000000000 0.000000e+00  0.000000e+00
# [7,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.01505735 -0.01505705  0.000000000 0.000000e+00  0.000000e+00
# [8,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -0.01505705  0.01703565  0.000000000 0.000000e+00  0.000000e+00
# [9,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.00000000  0.00000000  0.126293689 7.374380e-04 -1.248190e-01
# [10,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.00000000  0.00000000  0.000737438 3.744233e-04  9.216268e-07
# [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.00000000  0.00000000 -0.124819017 9.216268e-07  1.411815e-01
# 
# > testR$D_inv
# [[1]]
# NULL
# 
# [[2]]
# [,1]        [,2]        [,3]       [,4]       [,5]
# (Intercept)   1.284653   -3.137816   -37.08216  0.0000000  0.0000000
# year         -3.137816  106.671060  1372.13696  0.0000000  0.0000000
# I(year^2)   -37.082162 1372.136960 20472.52482  0.0000000  0.0000000
# (Intercept)   0.000000    0.000000     0.00000  0.2219492 -0.4376821
# year          0.000000    0.000000     0.00000 -0.4376821  9.4506461
# 
# [[3]]
# [,1]        [,2]        [,3]       [,4]         [,5]         [,6]       [,7]       [,8]
# (Intercept)   1.284653   -3.137816   -37.08216  0.0000000    0.0000000     0.000000  0.0000000  0.0000000
# year         -3.137816  106.671060  1372.13696  0.0000000    0.0000000     0.000000  0.0000000  0.0000000
# I(year^2)   -37.082162 1372.136960 20472.52482  0.0000000    0.0000000     0.000000  0.0000000  0.0000000
# (Intercept)   0.000000    0.000000     0.00000  1.8442427    0.1725615    -2.255727  0.0000000  0.0000000
# year          0.000000    0.000000     0.00000  0.1725615   41.2852252  -544.133635  0.0000000  0.0000000
# I(year^2)     0.000000    0.000000     0.00000 -2.2557275 -544.1336346 10327.658661  0.0000000  0.0000000
# (Intercept)   0.000000    0.000000     0.00000  0.0000000    0.0000000     0.000000  0.2219492 -0.4376821
# year          0.000000    0.000000     0.00000  0.0000000    0.0000000     0.000000 -0.4376821  9.4506461
# 
# [[4]]
# [,1]         [,2]         [,3]       [,4]       [,5]
# (Intercept)  1.8442427    0.1725615    -2.255727  0.0000000  0.0000000
# year         0.1725615   41.2852252  -544.133635  0.0000000  0.0000000
# I(year^2)   -2.2557275 -544.1336346 10327.658661  0.0000000  0.0000000
# (Intercept)  0.0000000    0.0000000     0.000000  0.2219492 -0.4376821
# year         0.0000000    0.0000000     0.000000 -0.4376821  9.4506461
# 
# 
# 

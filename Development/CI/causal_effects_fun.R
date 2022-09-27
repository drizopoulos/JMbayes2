if (FALSE) {
    object = jointFit1
    Data_Long = pt_all
    Data_Event = pt_all_CR
    t0 = 7
    Dt = 2
    IE_var = "salvage_indicator"
    IE_time = "salvage_time"
    B = 5L
    calculate_CI = TRUE
    cores = max(parallel::detectCores(logical = FALSE) - 1, 1)
    extra_objects = "dummy"
    seed = 1L

    Data <- get_data(pt_all, pt_all_CR, t0 = 8, Dt = 2, object = jointFit1,
                     IE_var = "salvage_indicator", IE_time = "salvage_time")
    Data_Long = Data$newdataL
    Data_Long2 = Data$newdataL2
    Data_Event = Data$newdataE
    Data_Event2 = Data$newdataE2


    Data_Long[Data_Long$pid == 1, ]
    Data_Long2[Data_Long2$pid == 1, ]

    Data_Event[Data_Event$pid == 1, c("pid", "start", "stop", "event",
                                      "salvage_indicator", "salvage_time", "CR")]

    Data_Event2[Data_Event2$pid == 1, c("pid", "start", "stop", "event",
                                      "salvage_indicator", "salvage_time", "CR")]

}


# a function to do the data management
get_data <- function (Data_Long, Data_Event, t0, Dt, object = NULL,
                      vars = NULL, IE_var = NULL, IE_time = NULL) {
    if (is.null(object) && is.null(vars)) {
        stop("one of the arguments 'object' or 'vars' must not be NULL.")
    }
    if (!is.null(object)) {
        if (!inherits(object, "jm")) {
            stop("this function works only with 'jm' model objects.")
        }
        if (object$model_info$type_censoring != "counting") {
            stop("this function only works using the counting process notation ",
                 "in the survival submodel.")
        }
        if (is.null(IE_var)) {
            stop("you also need to provide the 'IE_var'.")
        }
        if (is.null(IE_time)) {
            stop("you also need to provide the 'IE_time'.")
        }
        id_var <- object$model_info$var_names$idVar
        time_var <- object$model_info$var_names$time_var
        Time_var <- object$model_info$var_names$Time_var
        start_var <- Time_var[1L]
        Time_var <- Time_var[2L]
        event_var <- object$model_info$var_names$event_var
    } else {
        if (!is.character(vars)) {
            stop("'vars' must be a character vector.")
        }
        if (!(nams <- c('id_var', 'time_var', 'start_var', 'Time_var', 'event_var',
               'IE_var', 'IE_time')) %in% names(vars)) {
            stop("'vars' should be a named character vector with elements ",
                 " named: ", paste(nams, collapse = ", "))
        }
        id_var <- vars['id_var']
        time_var <- vars['time_var']
        start_var <- vars['start_var']
        Time_var <- vars['Time_var']
        event_var <- vars['event_var']
        # 'IE_var': character string with the name of the intermediate event variable.
        IE_var <- vars['IE_var']
        IE_time <- vars['IE_time']
    }
    ids_Long <- unique(Data_Long[[id_var]])
    ids_Event <- unique(Data_Event[[id_var]])
    if (!all(ids_Long %in% ids_Event) || !all(ids_Event %in% ids_Long)) {
        stop("the data.frames 'Data_Long' and 'Data_Event' ",
             "must contain the same subjects.")
    }
    # we find the subjects at risk at t0; for subjects who had an IE before
    # t0, we keep the measurement before the IE
    R_event <- Data_Event[Data_Event[[IE_var]] == 0, ]
    R_event <- R_event[R_event[[Time_var]] > t0, ]
    R_event[[id_var]] <- factor(R_event[[id_var]], unique(R_event[[id_var]]))
    # we set the stop time to t0
    R_event[[Time_var]] <- t0 + 1e-03
    # we set the event to zero
    R_event[[event_var]] <- 0
    # we create aof R_event with IE set to one
    R_event2 <- R_event
    R_event2[[IE_var]] <- 1
    R_event2[[IE_time]] <- t0

    # we keep the longitudinal measurements before t0 and before salvage for
    # patients who were at risk at t0
    check1 <- Data_Long[[time_var]] <= pmin(t0, Data_Long[[IE_time]])
    check2 <- Data_Long[[id_var]] %in% unique(R_event[[id_var]])
    R_long <- Data_Long[check1 & check2, ]
    R_long[[id_var]] <- factor(R_long[[id_var]])

    R_event <- R_event[R_event[[id_var]] %in% R_long[[id_var]], ]
    R_event[[id_var]] <- factor(R_event[[id_var]])

    R_event2 <- R_event2[R_event2[[id_var]] %in% R_long[[id_var]], ]
    R_event2[[id_var]] <- factor(R_event2[[id_var]])

    R_long2 <- R_long
    R_long2[[IE_var]] <- 1
    ###########################################################
    list(newdataL = R_long, newdataL2 = R_long2,
         newdataE = R_event, newdataE2 = R_event2)
}

causal_effects <- function (object, Data_Long, Data_Long2, Data_Event,
                            Data_Event2, t0, Dt, calculate_CI = FALSE, B = 50L,
                            cores = max(parallel::detectCores() - 1, 1),
                            extra_objects = NULL, seed = 1L) {
    # 'object': a joint model with survival submodel fitted with the
    # counting process notation and having a time-varying covariate indicating
    # an intermediate event.
    # 'Data_Long': a data.frame with the longitudinal measurements, containing
    # the intermediate event at zero
    # 'Data_Long2': a data.frame with the longitudinal measurements, containing
    # the intermediate event at one
    # 'Data_Event': a data.frame with the event times, containing the
    # intermediate event set at zero
    # 'Data_Event2': a data.frame with the event times, containing the
    # intermediate event set at one
    # 't0': numeric scalar denoting the time points at which to calculate
    # the causal effect.
    # 'Dt': numeric scalar of delta times, i.e., (t0, t0 + Dt)
    # 'calculate_CI': logical; should the function calculate the confidence
    # intervals of the effects.
    # 'B': numeric scalar denoting the number of Bootstrap samples.
    # 'cores': numeric scalar denoting the number of cores to be used.
    # 'extra_objects': character vector of objects to be passed to the workers.
    # 'seed': numeric scalar denoting the seed.
    ############################################################################
    # checks and extract variables
    if (!inherits(object, "jm")) {
        stop("this function works only with 'jm' model objects.")
    }
    if (object$model_info$type_censoring != "counting") {
        stop("this function only works using the counting process notation ",
             "in the survival submodel.")
    }
    if (!is.numeric(t0) || length(t0) > 1) {
        stop("'t0' must be a numeric scalar.")
    }
    if (!is.numeric(Dt) || length(Dt) > 1) {
        stop("'Dt' must be a numeric scalar.")
    }
    id_var <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    Time_var <- object$model_info$var_names$Time_var
    start_var <- Time_var[1L]
    Time_var <- Time_var[2L]
    event_var <- object$model_info$var_names$event_var
    vars <- c(id_var = id_var, time_var = time_var, start_var = start_var,
              Time_var = Time_var, event_var = event_var)
    all_vars_long <- unlist(lapply(object$model_info$terms$terms_FE, all.vars))
    all_vars_event <- all.vars(object$model_info$terms$terms_Surv)
    if (any(!c(all_vars_long, id_var) %in% names(Data_Long))) {
        stop("'Data_Long' should contain the variables: ",
             paste(unique(c(id_var, all_vars_long)), collapse = ", "))
    }
    if (any(!c(all_vars_event, id_var) %in% names(Data_Event))) {
        stop("'Data_Event' should contain the variables: ",
             paste(unique(c(id_var, all_vars_event)), collapse = ", "))
    }
    if (any(!c(all_vars_event, id_var) %in% names(Data_Event2))) {
        stop("'Data_Event2' should contain the variables: ",
             paste(unique(c(id_var, all_vars_event)), collapse = ", "))
    }
    ids_Long <- unique(Data_Long[[id_var]])
    ids_Long2 <- unique(Data_Long2[[id_var]])
    ids_Event <- unique(Data_Event[[id_var]])
    ids_Event2 <- unique(Data_Event2[[id_var]])
    if (!all(ids_Long %in% ids_Event) || !all(ids_Long %in% ids_Event2)
        || !all(ids_Event %in% ids_Long) || !all(ids_Event2 %in% ids_Long)
        || !all(ids_Event %in% ids_Event2) || !all(ids_Event2 %in% ids_Event)
        || !all(ids_Long2 %in% ids_Event) || !all(ids_Event2 %in% ids_Long2)) {
        stop("the data.frames 'Data_Long', 'Data_Long2', 'Data_Event' and 'Data_Event2' ",
             "must contain the same subjects.")
    }
    # a function to calculate the causal effects (per stratum) in the presence
    # and absence of the IE
    get_effect <- function (object, Data_Long, Data_Long2, Data_Event,
                            Data_Event2, t0, Dt, vars) {
        # we calculate the CIFs without IE
        newdata_withoutIE <- list(newdataL = Data_Long, newdataE = Data_Event)

        CIF_withoutIE <- predict(object, newdata = newdata_withoutIE,
                                 process = "event", times = t0 + Dt,
                                 return_mcmc = TRUE)

        # we calculate the CIFs with IE
        newdata_withIE <- list(newdataL = Data_Long2, newdataE = Data_Event2)
        CIF_withIE <- predict(object, newdata = newdata_withoutIE,
                              newdata2 = newdata_withIE,
                              process = "event", times = t0 + Dt,
                              return_mcmc = TRUE)

        # the marginal effect is the mean over the conditional effects
        ind <- CIF_withIE$times > t0
        strata <- CIF_withIE[["_strata"]][ind]
        cif1 <- CIF_withIE$pred[ind]
        cif2 <- CIF_withoutIE$pred[ind]
        # drop failed iters
        cif1[cif1 > 1 | cif1 < 0] <- as.numeric(NA)
        cif2[cif2 > 1 | cif2 < 0] <- as.numeric(NA)
        effect <- tapply(cif1 - cif2, strata, mean, na.rm = TRUE)
        mcmc1 <- CIF_withIE$mcmc[ind, ]
        mcmc2 <- CIF_withoutIE$mcmc[ind, ]
        mcmc1[mcmc1 > 1 | mcmc1 < 0] <- as.numeric(NA)
        mcmc2[mcmc2 > 1 | mcmc2 < 0] <- as.numeric(NA)
        vv <- rowsum(mcmc1 - mcmc2, strata, reorder = FALSE, na.rm = TRUE) /
            sum(strata == 1)
        attr(effect, "var") <- matrixStats::rowVars(vv)
        effect
    }
    # a function to create a non-parametric Bootstrap sample
    make_bootSample <- function (Data_Long, Data_Long2, Data_Event, Data_Event2,
                                 id_var) {
        ids <- Data_Long[[id_var]]
        unq_ids <- unique(ids)
        ids <- factor(ids, levels = unq_ids)
        new_ids <- sample(unq_ids, replace = TRUE)
        new_Data_Long <- new_Data_Long2 <- new_Data_Event <- new_Data_Event2 <-
            vector("list", length(unq_ids))
        for (i in seq_along(unq_ids)) {
            keep <- Data_Long[[id_var]] == new_ids[i]
            dataL_i <- Data_Long[keep, ]
            dataL_i[[id_var]] <- i
            new_Data_Long[[i]] <- dataL_i
            ##
            keep <- Data_Long2[[id_var]] == new_ids[i]
            dataL2_i <- Data_Long2[keep, ]
            dataL2_i[[id_var]] <- i
            new_Data_Long2[[i]] <- dataL2_i
            ##
            keep <- Data_Event[[id_var]] == new_ids[i]
            dataE_i <- Data_Event[keep, ]
            dataE_i[[id_var]] <- i
            new_Data_Event[[i]] <- dataE_i
            ##
            keep <- Data_Event2[[id_var]] == new_ids[i]
            dataE2_i <- Data_Event2[keep, ]
            dataE2_i[[id_var]] <- i
            new_Data_Event2[[i]] <- dataE2_i

        }
        list(Data_Long = do.call("rbind", new_Data_Long),
             Data_Long2 = do.call("rbind", new_Data_Long2),
             Data_Event = do.call("rbind", new_Data_Event),
             Data_Event2 = do.call("rbind", new_Data_Event2))
    }
    #######################################################
    # causal effects in the original data
    effects <- get_effect(object, Data_Long, Data_Long2, Data_Event,
                          Data_Event2, t0, Dt, vars)
    if (calculate_CI) {
        # run Bootstrap in parallel
        if (!exists(".Random.seed", envir = .GlobalEnv))
            runif(1L)
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
        samples <- split(seq_len(B), rep(seq_len(cores), each = ceiling(B / cores),
                                         length.out = B))
        boot_parallel <- function (samples, object, Data_Long, Data_Long2,
                                   Data_Event, Data_Event2, t0, Dt, vars) {
            str <- object$model_data$strata
            out <- matrix(0.0, length(samples), length(unique(str)))
            for (b in seq_along(samples)) {
                boot <- make_bootSample(Data_Long, Data_Long2,
                                        Data_Event, Data_Event2, vars['id_var'])
                meffects <- get_effect(object, boot$Data_Long, boot$Data_Long2,
                                       boot$Data_Event,
                                       boot$Data_Event2, t0, Dt, vars)
                out[b, ] <- meffects
            }
            out
        }
        if (cores > 1L) {
            cl <- parallel::makeCluster(cores)
            parallel::clusterExport(cl, c("make_bootSample", "get_effect", extra_objects),
                                    envir = environment())
            parallel::clusterEvalQ(cl = cl, library("JMbayes2"))
            parallel::clusterSetRNGStream(cl = cl, iseed = seed)
            out <- parallel::parLapply(cl, samples, boot_parallel,
                                       object = object, Data_Long = Data_Long,
                                       Data_Long2 = Data_Long2,
                                       Data_Event = Data_Event,
                                       Data_Event2 = Data_Event2,
                                       t0 = t0, Dt = Dt, vars = vars)
            parallel::stopCluster(cl)
        } else {
            out <- lapply(samples, boot_parallel,
                          object = object, Data_Long = Data_Long,
                          Data_Long2 = Data_Long2,
                          Data_Event = Data_Event, Data_Event2 = Data_Event2,
                          t0 = t0, Dt = Dt, vars = vars)
        }
        out <- do.call('rbind', out)
        var_effects <- matrixStats::colVars(out) + attr(effects, "var")
        names(var_effects) <- names(effects)
    }
    attr(effects, "var") <- NULL
    list("effects" = effects,
         "CIs" = if (calculate_CI)
             cbind(low = effects - 1.96 * sqrt(var_effects),
                   upp = effects + 1.96 * sqrt(var_effects)))
}

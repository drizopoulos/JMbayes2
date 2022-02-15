object = jointFit
Data_Long = pbc2
Data_Event = pbc2_CR
t0 = 3
Dt = 2
IE_var = "IE"
B = 50L
cores = max(parallel::detectCores(logical = FALSE) - 1, 1)
extra_objects = "dummy"
seed = 1L

causal_effects <- function (object, Data_Long, Data_Event, t0, Dt, IE_var,
                            B = 50L, cores = max(parallel::detectCores() - 1, 1),
                            extra_objects = NULL, seed = 1L) {
    # 'object': a joint model with survival submodel fitted with the
    # counting process notation and having a time-varying covariate indicating
    # an intermediate event.
    # 'Data_Long': a data.frame with the longitudinal measurements, containing
    # the intermediate event
    # 'Data_Event': a data.frame with the event times, containing the
    # intermediate event
    # 't0': numeric scalar denoting the time points at which to calculate
    # the causal effect.
    # 'Dt': numeric scalar of delta times, i.e., (t0, t0 + Dt)
    # 'IE_var': character string with the name of the intermediate event variable.
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
    id_var <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    Time_var <- object$model_info$var_names$Time_var
    start_var <- Time_var[1L]
    Time_var <- Time_var[2L]
    event_var <- object$model_info$var_names$event_var
    vars <- c(id_var = id_var, time_var = time_var, start_var = start_var,
              Time_var = Time_var, event_var = event_var, IE_var = IE_var)
    all_vars_long <- unlist(lapply(object$model_info$terms$terms_FE, all.vars))
    all_vars_event <- all.vars(object$model_info$terms$terms_Surv)
    if (any(!c(all_vars_long, id_var, IE_var) %in% names(Data_Long))) {
        stop("'Data_Long' should contain the variables: ",
             paste(unique(c(id_var, all_vars_long, IE_var)), collapse = ", "))
    }
    if (any(!c(all_vars_event, id_var, IE_var) %in% names(Data_Event))) {
        stop("'Data_Event' should contain the variables: ",
             paste(unique(c(id_var, all_vars_event, IE_var)), collapse = ", "))
    }
    # a function to do the data management
    get_data <- function (Data_Long, Data_Event, t0, Dt, vars) {
        id_var <- vars['id_var']
        time_var <- vars['time_var']
        start_var <- vars['start_var']
        Time_var <- vars['Time_var']
        event_var <- vars['event_var']
        IE_var <- vars['IE_var']
        # we find the subjects at risk at t0
        keep <- function (x) rep(max(x) > t0, length(x))
        at_risk <- ave(Data_Event[[Time_var]], Data_Event[[id_var]], FUN = keep)
        R_event <- Data_Event[as.logical(at_risk), ]
        R_event[[id_var]] <- factor(R_event[[id_var]], unique(R_event[[id_var]]))
        # then we want to keep only the subjects who did not have an IE
        # up to time t0
        spl <- split(R_event, R_event[[id_var]])
        keep <- function (d) {
            start_times <- d[[start_var]]
            IEs <- d[[IE_var]]
            v <- if (all(IEs == 0)) TRUE else max(start_times[IEs == 1]) > t0
            rep(v, nrow(d))
        }
        R_event <- R_event[unlist(lapply(spl, keep)), ]
        R_event[[id_var]] <- factor(R_event[[id_var]], unique(R_event[[id_var]]))
        # there are some subjects who had an IE after t0; we need
        # to remove these lines from R_event
        R_event <- R_event[R_event[[IE_var]] == 0, ]
        # we set the stop time to t0
        R_event[[Time_var]] <- t0
        # we set the event to zero
        R_event[[event_var]] <- 0
        # we create aof R_event with IE set to one
        R_event2 <- R_event
        R_event2[[IE_var]] <- 1

        # we keep the longitudinal measurements before t0 for patients who
        # were at risk at t0
        check1 <- Data_Long[[time_var]] <= t0
        check2 <- Data_Long[[id_var]] %in% unique(R_event[[id_var]])
        R_long <- Data_Long[check1 & check2, ]
        R_long[[id_var]] <- factor(R_long[[id_var]])
        ###########################################################
        list(newdataL = R_long, newdataE = R_event, newdataE2 = R_event2)
    }
    # a function to calculate the causal effect in the presence and absence of
    # the IE
    get_effect <- function (object, Data_Long, Data_Event, t0, Dt, vars) {
        Data_i <- get_data(Data_Long, Data_Event, t0, Dt, vars)
        # we calculate the CIFs without IE
        newdata_withoutIE <- list(newdataL = Data_i[["newdataL"]],
                                  newdataE = Data_i[["newdataE"]])

        CIF_withoutIE <- predict(object, newdata = newdata_withoutIE,
                                 process = "event", times = t0 + Dt,
                                 return_mcmc = TRUE)

        # we calculate the CIFs with IE
        newdata_withIE <- list(newdataL = Data_i[["newdataL"]],
                               newdataE = Data_i[["newdataE2"]])
        CIF_withIE <- predict(object, newdata = newdata_withoutIE,
                              newdata2 = newdata_withIE,
                              process = "event", times = t0 + Dt,
                              return_mcmc = TRUE)

        # the marginal effect is the mean over the conditional effects
        mean(CIF_withIE$pred[CIF_withIE$times > t0] -
                 CIF_withoutIE$pred[CIF_withoutIE$times > t0])
    }
    # a function to create a Bootstrap sample
    make_bootSample <- function (Data_Long, Data_Event, id_var) {
        ids <- Data_Long[[id_var]]
        unq_ids <- unique(ids)
        ids <- factor(ids, levels = unq_ids)
        new_ids <- sample(unq_ids, replace = TRUE)
        new_Data_Long <- new_Data_Event <- vector("list", length(unq_ids))
        for (i in seq_along(unq_ids)) {
            keep <- Data_Long[[id_var]] == new_ids[i]
            dataL_i <- Data_Long[keep, ]
            dataL_i[[id_var]] <- i
            new_Data_Long[[i]] <- dataL_i
            ##
            keep <- Data_Event[[id_var]] == new_ids[i]
            dataE_i <- Data_Event[keep, ]
            dataE_i[[id_var]] <- i
            new_Data_Event[[i]] <- dataE_i
        }
        list(Data_Long = do.call("rbind", new_Data_Long),
             Data_Event = do.call("rbind", new_Data_Event))
    }
    #######################################################
    # causal effect in the original data
    effect <- get_effect(object, Data_Long, Data_Event, t0, Dt, vars)
    # run Bootstrap
    if (!exists(".Random.seed", envir = .GlobalEnv))
        runif(1L)
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", RNGstate, envir = .GlobalEnv))
    samples <- split(seq_len(B), rep(seq_len(cores), each = ceiling(B / cores),
                                     length.out = B))
    boot_parallel <- function (samples, object, Data_Long, Data_Event, t0,
                               Dt, vars) {
        out <- numeric(length(samples))
        for (b in seq_along(samples)) {
            boot <- make_bootSample(Data_Long, Data_Event, vars['id_var'])
            meffect <- get_effect(object, boot$Data_Long, boot$Data_Event,
                                  t0, Dt, vars)
            out[b] <- meffect
        }
        out
    }
    cl <- parallel::makeCluster(cores)
    parallel::clusterExport(cl, c("make_bootSample", "get_effect", "get_data",
                                  extra_objects), envir = environment())
    parallel::clusterEvalQ(cl = cl, library("JMbayes2"))
    parallel::clusterSetRNGStream(cl = cl, iseed = seed)
    out <- parallel::parLapply(cl, samples, boot_parallel,
                               object = object, Data_Long = Data_Long,
                               Data_Event = Data_Event, t0 = t0,
                               Dt = Dt, vars = vars)
    parallel::stopCluster(cl)
    out <- unlist(out)
    list("effect" = effect, "var_effect" = var(out))
}

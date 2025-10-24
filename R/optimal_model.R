#if (FALSE) {
#    library("JMbayes2")
#    pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
#    pbc2$status2 <- as.numeric(pbc2$status != 'alive')
#    CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
#    fm1 <- lme(log(serBilir) ~ ns(year, 3) * sex, data = pbc2,
#               random = list(id = pdDiag(~ ns(year, 3))))
#    fm2 <- lme(log(serBilir) ~ poly(year, 2) * sex, data = pbc2,
#               random = list(id = pdDiag(~ poly(year, 2))))
#    jointFit1 <- jm(CoxFit, list(fm1), time_var = "year")
#    jointFit2 <- jm(CoxFit, list(fm2), time_var = "year")
#    #####
#    models = list(jointFit1, jointFit2)
#    t0 = 8.5
#    newdata = pbc2[ave(pbc2$year, pbc2$id, FUN = max) > t0, ]
#}

opt_model <- function (models, newdata, t0, parallel = "snow", cores = 1L) {
    if (!all(sapply(models, function (obj) inherits(obj, "jm")))) {
        stop("all objects in 'models' must inherit from class 'jm'.\n")
    }
    object <- models[[1L]]
    id_var <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    Time_var <- object$model_info$var_names$Time_var
    TTime_var <- if (length(Time_var) > 1) Time_var[2L] else Time_var
    event_var <- object$model_info$var_names$event_var
    if (!is.data.frame(newdata)) {
        stop("'newdata' must be a data.frame.\n")
    }
    if (!all((c(id_var, time_var) %in% names(newdata)))) {
        stop("'newdata' must contain the variables '", id_var, "' and '",
             time_var, "'.\n")
    }
    ND <- newdata
    ND[[id_var]] <- match(ND[[id_var]], unique(ND[[id_var]]))
    ND_before <- ND[ND[[time_var]] <= t0, ]
    ND_before[[TTime_var]] <- t0
    ND_before[[event_var]] <- 0
    ND_after <- ND[ND[[time_var]] > t0, ]
    data_after_t0 <- nrow(ND_after)
    ###
    MISE_perid <- function (obs, reps, id, times) {
        trapezoid_rule <- function (f, x) {
            sum(0.5 * diff(x) * (f[-length(x)] + f[-1L]))
        }
        loess.smooth2 <- function (x, y) {
            loess.smooth(x, y, degree = 2, span = 0.75,
                         family = "gaussian", evaluation = 200)
        }
        smooth <- function (x, y) {
            n <- length(x)
            if (n > 5) {
                loess.smooth2(x, y)
            } else if (n > 1 && n <= 5) {
                approx(x, y)
            } else {
                list(x = NA_real_, y = NA_real_)
            }
        }
        gof_fun <- function (y, times, id, type) {
            if (type == "variogram") {
                ls <- loess.smooth2(times, y)
                ind <- findInterval(times, ls$x)
                rr <- y - ls$y[ind]
                variogram(rr, times, id)[[1L]]
            } else if (type == "variance-function") {
                ls <- loess.smooth2(times, y)
                ind <- findInterval(times, ls$x)
                rr <- y - ls$y[ind]
                sigma <- sqrt(sum(rr * rr) / (length(rr) - 5.35))
                rr <- sqrt(abs(rr / sigma))
                cbind(times, rr)
            } else {
                cbind(times, y)
            }
        }
        Obs_ave <- gof_fun(obs, times, id, "average")
        Obs_vario <- gof_fun(obs, times, id, "variogram")
        F_obs_ave <- mapply(smooth, y = split(Obs_ave[, 2L], id),
                            x = split(Obs_ave[, 1L], id), SIMPLIFY = FALSE)
        ni <- tapply(id, id, length)
        id_long <- rep(unique(id), sapply(ni, function (n) ncol(combn(n, 2))))
        F_obs_vario <- mapply(smooth, y = split(Obs_vario[, 2L], id_long),
                              x = split(Obs_vario[, 1L], id_long),
                              SIMPLIFY = FALSE)
        mise <- function (obs, rep) {
            trapezoid_rule((obs$y - rep$y)^2, obs$x)
        }
        n <- length(unique(id))
        M <- ncol(reps)
        MISE <- matrix(0.0, n, M)
        for (m in seq_len(M)) {
            reps_ave <- gof_fun(reps[, m], times, id, "average")
            reps_vario <- gof_fun(reps[, m], times, id, "variogram")
            F_reps_ave <-
                mapply(smooth, y = split(reps_ave[, 2L], id),
                       x = split(reps_ave[, 1L], id), SIMPLIFY = FALSE)
            F_reps_vario <-
                mapply(smooth, y = split(reps_vario[, 2L], id_long),
                       x = split(reps_vario[, 1L], id_long), SIMPLIFY = FALSE)
            mise_ave <- mapply(mise, obs = F_obs_ave, rep = F_reps_ave)
            mise_vario <- mapply(mise, obs = F_obs_vario, rep = F_reps_vario)
            MISE[, m] <- c(scale(mise_ave)) + c(scale(mise_vario))
        }
        rowMeans(MISE)
    }
    MISE_model <- function (object) {
        prs <- predict(object, newdata = ND_before, newdata2 = ND_after,
                       return_params_mcmc = TRUE)
        sims <- simulate(object, nsim = 200L, newdata = ND_before,
                          include_outcome = TRUE, random_effects = "mcmc",
                          params_mcmc = prs$newdata$params_mcmc)
        n_outcomes <- length(sims[["outcome"]])
        id <- attr(sims$outcome[[1]], "id")
        MISEs <- matrix(0, length(unique(id)), n_outcomes)
        colnames(MISEs) <- names(sims)[seq_len(n_outcomes)]
        for (m in seq_len(n_outcomes)) {
            outcome <- sims$outcome[[m]]
            MISEs[, m] <-
                MISE_perid(outcome, sims[[m]], attr(outcome, "id"),
                           attr(outcome, "times"))
        }
        list(MISEs = MISEs, Preds = prs)
    }
    ###
    if (cores > 1L && length(models) > 1L) {
        have_mc <- have_snow <- FALSE
        if (parallel == "multicore") {
            have_mc <- .Platform$OS.type != "windows"
        } else if (parallel == "snow") {
            have_snow <- TRUE
        }
        if (!have_mc && !have_snow) cores <- 1L
        loadNamespace("parallel")
        if (have_mc) {
            out <- parallel::mclapply(models, MISE_model, mc.cores = cores)
        } else {
            cl <- parallel::makePSOCKcluster(rep("localhost", cores))
            invisible(parallel::clusterEvalQ(cl, library("JMbayes2")))
            out <- parallel::parLapply(cl, models, MISE_model)
            parallel::stopCluster(cl)
        }
    } else {
        out <- lapply(models, MISE_model)
    }
    out[] <- lapply(out, function (x, ids) {rownames(x[[1L]]) <- ids; x},
                    ids = unique(newdata[[id_var]]))
    out
}






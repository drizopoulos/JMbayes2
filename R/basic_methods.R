traceplot <- function (object, ...) UseMethod("traceplot")

traceplot.jm <- function (object,
                          parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                                   "tau_bs_gammas", "gammas", "alphas"),
                          ...) {
    parm <- match.arg(parm)
    if (parm == "all") {
        parms <- c("betas", "sigmas", "D", "bs_gammas", "tau_bs_gammas",
                   "gammas", "alphas")
        for (i in seq_along(parms)) {
            parms_i <- parms[[i]]
            x <- object$mcmc[[parms_i]]
            if (!is.null(x)) coda::traceplot(x, ...)
        }
    } else {
        coda::traceplot(object$mcmc[[parm]], ...)
    }
    invisible()
}

gelman_diag <- function (object, ...) UseMethod("gelman_diag")

gelman_diag.jm <- function (object,
                          parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                                   "tau_bs_gammas", "gammas", "alphas"),
                          ...) {
    parm <- match.arg(parm)
    if (parm == "all") {
        parms <- c("betas", "sigmas", "D", "bs_gammas", "tau_bs_gammas",
                   "gammas", "alphas")
        out <- vector("list", length(parms))
        names(out) <- parms
        for (i in seq_along(parms)) {
            parms_i <- tail(grep(parms[[i]], names(object$mcmc)), 1L)
            x <- object$mcmc[[parms_i]]
            if (!is.null(x)) out[[i]] <- coda::gelman.diag(x, ...)
        }
        out[!sapply(out, is.null)]
    } else {
        parm <- tail(grep(parm, names(object$mcmc)), 1L)
        coda::gelman.diag(object$mcmc[[parm]], ...)
    }
}

densityplot <- function (object, ...) UseMethod("densityplot")

densityplot.jm <- function (object,
                          parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                                   "tau_bs_gammas", "gammas", "alphas"),
                          ...) {
    parm <- match.arg(parm)
    if (parm == "all") {
        parms <- c("betas", "sigmas", "D", "bs_gammas", "tau_bs_gammas",
                   "gammas", "alphas")
        for (i in seq_along(parms)) {
            parms_i <- parms[[i]]
            x <- object$mcmc[[parms_i]]
            if (!is.null(x)) coda::densityplot(x, ...)
        }
    } else {
        coda::densityplot(object$mcmc[[parm]], ...)
    }
    invisible()
}

summary.jm <- function (object, ...) {
    families <- object$model_info$families
    n_outcomes <- length(families)
    respVars <- object$model_info$var_names$respVars_form
    N <- sapply(object$model_data$X, nrow)
    descrpt <- data.frame(` ` = N, row.names = respVars, check.rows = FALSE,
                          check.names = FALSE)
    nams_D <- unlist(lapply(object$model_data$Z, colnames))
    D <- lowertri2mat(object$statistics$Mean$D, nams_D)
    out <- list(n = object$model_data$n, descrpt = descrpt, D = D,
                families = families, respVars = respVars,
                events = object$model_data$delta,
                control = object$control, time = object$running_time,
                call = object$call)
    tab_f <- function(name) {
        out <- data.frame(Mean = object$statistics$Mean[[name]],
                          StDev = object$statistics$SD[[name]],
                          `2.5%` = object$statistics$CI_low[[name]],
                          `97.5%` = object$statistics$CI_low[[name]],
                          P = object$statistics$P[[name]],
                          row.names = names(object$statistics$P[[name]]),
                          check.names = FALSE)
        Rhat <- object$statistics$Rhat[[name]][, 1L]
        if (!is.null(Rhat))
            out$Rhat <- Rhat
        out
    }
    fam_names <- sapply(families, "[[", "family")
    has_sigma <- c("gaussian", "Student-t", "beta", "Gamma",
                   "negative binomial", "beta binomial")
    tab_sigmas <- tab_f("sigmas")
    for (i in seq_len(n_outcomes)) {
        nam_outcome <- paste0("Outcome", i)
        out[[nam_outcome]] <- tab_f(paste0("betas", i))
        if (fam_names[i] %in% has_sigma) {
            k <- nrow(out[[nam_outcome]])
            out[[nam_outcome]] <- rbind(out[[nam_outcome]], tab_sigmas[i, ])
            row.names(out[[nam_outcome]])[k + 1] <- "sigma"
        }
    }
    out$Survival <- do.call(rbind, list(tab_f("gammas"), tab_f("alphas")))
    class(out) <- "summary.jm"
    out
}

print.summary.jm <- function (x, digits = max(4, getOption("digits") - 4), ...) {
    cat("\nCall:\n", printCall(x$call), "\n\n", sep = "")
    cat("Data Descriptives:")
    cat("\nNumber of Groups: ", x$n, "\t\tNumber of events: ",
        sum(x$event == 1), " (", round(100 * mean(x$event == 1), 1),
        "%)", sep = "")
    cat("\nNumber of Observations:")
    obs <- x$descrpt
    for (i in 1:nrow(obs)) {
        cat("\n  ", row.names(obs)[i], ": ", obs[[1]][i],
            sep = "")
    }
    cat("\n")
    if (!is.null(x$DIC)) {
        model.sum <- data.frame(DIC = x$DIC, pD = x$pD, row.names = "")
        print(model.sum)
    }
    cat("\nRandom-effects covariance matrix:\n")
    D <- x$D
    ncz <- nrow(D)
    sds <- sqrt(diag(D))
    if (ncz > 1) {
        corrs <- cov2cor(D)
        corrs[upper.tri(corrs, TRUE)] <- 0
        mat <- round(cbind(sds, corrs[, -ncz]), digits)
        mat <- rbind(mat)
        mat <- apply(mat, 2L, sprintf, fmt = "%.4f")
        mat[mat == mat[1, 2]] <- ""
        mat[1, -1] <- sprintf("%06s", abbreviate(colnames(D)[-ncz], 6))
        colnames(mat) <- rep("", ncol(mat))
        mat <- rbind(c("StdDev", "  Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL),
                     mat)
        rownames(mat) <- c("", abbreviate(c(dimnames(D)[[1]]), 6))
    } else {
        mat <- cbind(StdDev = sprintf(sds, fmt = "%.4f"))
        rownames(mat) <- rownames(D)
    }
    print(noquote(mat), digits = digits)
    cat("\nSurvival Outcome:\n")
    print(round(x[["Survival"]], digits))
    n_outcomes <- length(x$families)
    for (i in seq_len(n_outcomes)) {
        cat("\nLongitudinal Outcome: ", x$respVars[i],
            " (family = ", x$families[[i]][["family"]],
            ", link = ", x$families[[i]][["link"]],
            ")", "\n", sep = "")
        xx <- round(x[[paste0("Outcome", i)]], digits)
        rnams <- row.names(xx)
        if (any(offend <- nchar(rnams) > 20))
            row.names(xx)[offend] <- abbreviate(rnams[offend])
        print(xx)
    }
    cat("\nMCMC summary:\n")
    tt <- x$time[3L] / 60
    cat("chains:", x$control$n_chains,
        "\niterations per chain:", x$control$n_iter,
        "\nburn-in per chain:", x$control$n_burnin,
        "\ntime:", if (tt > 60)
            round(tt/60, 1)
        else round(tt, 1), if (tt > 60)
            "hours"
        else "min")
    cat("\n")
    invisible(x)
}

print.jm <- function (x, digits = max(4, getOption("digits") - 4), ...) {
    xx <- summary(x)
    cat("\nCall:\n", printCall(xx$call), "\n", sep = "")
    cat("\nRandom-effects covariance matrix:\n")
    D <- xx$D
    ncz <- nrow(D)
    diag.D <- ncz != ncol(D)
    sds <- if (diag.D) sqrt(D) else sqrt(diag(D))
    if (ncz > 1) {
        if (diag.D) {
            dat <- as.data.frame(round(rbind(sds), digits))
            names(dat) <- "StdDev"
        } else {
            corrs <- cov2cor(D)
            corrs[upper.tri(corrs, TRUE)] <- 0
            mat <- round(cbind(sds, corrs[, -ncz]), digits)
            mat <- rbind(mat)
            mat <- apply(mat, 2, sprintf, fmt = "% .4f")
            mat[mat == mat[1, 2]] <- ""
            mat[1, -1] <- abbreviate(colnames(D)[-ncz], 6)
            colnames(mat) <- c(colnames(mat)[1], rep("",
                                                     ncz - 1))
            dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
            names(dat) <- c("StdDev", "Corr", if (ncz >
                                                  2) rep(" ", ncz - 2) else NULL)
            row.names(dat) <- abbreviate(c(dimnames(D)[[1]]))
        }
    } else {
        dat <- data.frame(StdDev = c(sds, x$sigma), row.names = if (!is.null(x$sigma))
            c(rownames(D), "Residual")
            else rownames(D), check.rows = FALSE, check.names = FALSE)
    }
    print(dat, digits = digits)
    cat("\nSurvival Outcome:\n")
    print(round(xx[["Survival"]], digits))
    n_outcomes <- length(xx$families)
    for (i in seq_len(n_outcomes)) {
        cat("\nLongitudinal Outcome: ", xx$respVars[i],
            " (family = ", xx$families[[i]][["family"]],
            ", link = ", xx$families[[i]][["link"]],
            ")", "\n", sep = "")
        yy <- round(xx[[paste0("Outcome", i)]], digits)
        rnams <- row.names(yy)
        if (any(offend <- nchar(rnams) > 20))
            row.names(yy)[offend] <- abbreviate(rnams[offend])
        print(yy)
    }
    cat("\n")
    invisible(x)
}

fixef.jm <- function(object, outcome = 1, ...) {
    if (!is.numeric(outcome) || outcome < 0) {
        stop("'outcome' should be a positive integer.")
    }
    outcome <- round(outcome)
    if (outcome > length(object$model_data$y)) {
        Means <- object$statistics$Mean
        ind_betas <- grep("betas", names(Means), fixed = TRUE)
        Means <- Means[ind_betas]
        names(Means) <- object$model_info$var_names$respVars_form
        Means
    } else {
        object$statistics$Mean[[paste0("betas", outcome)]]
    }
}

ranef.jm <- function(object, outcome = Inf, ...) {
    if (!is.numeric(outcome) || outcome < 0) {
        stop("'outcome' should be a positive integer.")
    }
    outcome <- round(outcome)
    if (outcome > length(object$model_data$y)) {
        object$statistics$RE
    } else {
        Means <- object$statistics$Mean$RE
        ind <- 1:2 # to be fixed.
        Means[, ind]
    }
}

terms.MixMod <- function (x, type = c("fixed", "random", "zi_fixed", "zi_random"), ...) {
    type <- match.arg(type)
    switch(type, "fixed" = x$Terms$termsX, "random" = x$Terms$termsZ,
           "zi_fixed" = x$Terms$termsX_zi, "zi_random" = x$Terms$termsZ_zi)
}

model.frame.MixMod <- function (formula, type = c("fixed", "random", "zi_fixed",
                                                  "zi_random"), ...) {
    type <- match.arg(type)
    switch(type, "fixed" = formula$model_frames$mfX, "random" = formula$model_frames$mfZ,
           "zi_fixed" = formula$model_frames$mfX_zi,
           "zi_random" = formula$model_frames$mfZ_zi)
}

model.matrix.MixMod <- function (object, type = c("fixed", "random", "zi_fixed", "zi_random"), ...) {
    type <- match.arg(type)
    switch(type,
           "fixed" = model.matrix(object$Terms$termsX, object$model_frames$mfX),
           "random" = {
               id <- object$id[[1]]
               id <- match(id, unique(id))
               Z <- mapply(constructor_Z, object$Terms$termsZ, object$model_frames$mfZ,
                           MoreArgs = list(id = id), SIMPLIFY = FALSE)
               do.call("cbind", Z)
           },
           "zi_fixed" = model.matrix(object$Terms$termsX_zi, object$model_frames$mfX_zi),
           "zi_random" = {
               id <- object$id[[1]]
               id <- match(id, unique(id))
               Z <- mapply(constructor_Z, object$Terms$termsZ_zi, object$model_frames$mfZ_zi,
                           MoreArgs = list(id = id), SIMPLIFY = FALSE)
               do.call("cbind", Z)
           }
    )
}

formula.MixMod <- function (x, type = c("fixed", "random", "zi_fixed", "zi_random"), ...) {
    type <- match.arg(type)
    switch(type, "fixed" = eval(x$call$fixed), "random" = eval(x$call$random),
           "zi_fixed" = eval(x$call$zi_fixed), "zi_random" = eval(x$call$zi_random))
}

family.MixMod <- function (object, ...) {
    object$family
}

nobs.MixMod <- function (object, level = 1,...) {
    if (level == 0) {
        length(unique(object$id[[1]]))
    } else {
        length(object$id[[1]])
    }
}






traceplot <- function (object, ...) UseMethod("traceplot")

traceplot.jm <- function (object,
                          parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                                   "tau_bs_gammas", "gammas", "alphas"),
                          ...) {
    parm <- match.arg(parm)
    if (parm == "all") {
        nams_parms <- c("betas", "sigmas", "D", "bs_gammas", "tau_bs_gammas",
                        "gammas", "alphas")
        nams_mcmc <- names(object$mcmc)
        ind <- unlist(sapply(paste0("^", nams_parms), grep, nams_mcmc))
        nams_mcmc <- nams_mcmc[ind]
        for (i in seq_along(nams_mcmc)) {
            parms_i <- nams_mcmc[[i]]
            x <- object$mcmc[[parms_i]]
            if (!is.null(x)) coda::traceplot(x, ...)
        }
    } else {
        parm <- grep(paste0("^", parm), names(object$mcmc))
        if (length(parm) > 1) {
            for (l in parm) coda::traceplot(object$mcmc[[l]], ...)
        } else {
            coda::traceplot(object$mcmc[[parm]], ...)
        }
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
        nams_parms <- c("betas", "sigmas", "D", "bs_gammas", "tau_bs_gammas",
                        "gammas", "alphas")
        nams_mcmc <- names(object$mcmc)
        ind <- unlist(sapply(paste0("^", nams_parms), grep, nams_mcmc), use.names= FALSE)
        nams_mcmc <- nams_mcmc[ind]
        out <- vector("list", length(nams_mcmc))
        names(out) <- nams_mcmc
        for (i in seq_along(out)) {
            parms_i <- nams_mcmc[[i]]
            x <- object$mcmc[[parms_i]]
            if (!is.null(x)) out[[i]] <- coda::gelman.diag(x, ...)
        }
        out[!sapply(out, is.null)]
    } else {
        parm <- grep(paste0("^", parm), names(object$mcmc))
        if (length(parm) > 1) {
            out <- lapply(parm, function (l)
                coda::gelman.diag(object$mcmc[[l]], ...))
            names(out) <- object$model_info$var_names$respVars_form
            out
        } else {
            coda::gelman.diag(object$mcmc[[parm]], ...)
        }
    }
}

densplot <- function (object, ...) UseMethod("densplot")

densplot.jm <- function (object,
                          parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                                   "tau_bs_gammas", "gammas", "alphas"),
                          ...) {
    parm <- match.arg(parm)
    if (parm == "all") {
        nams_parms <- c("betas", "sigmas", "D", "bs_gammas", "tau_bs_gammas",
                        "gammas", "alphas")
        nams_mcmc <- names(object$mcmc)
        ind <- unlist(sapply(paste0("^", nams_parms), grep, nams_mcmc))
        nams_mcmc <- nams_mcmc[ind]
        for (i in seq_along(nams_mcmc)) {
            parms_i <- nams_mcmc[[i]]
            x <- object$mcmc[[parms_i]]
            if (!is.null(x)) coda::densplot(x, ...)
        }
    } else {
        parm <- grep(paste0("^", parm), names(object$mcmc))
        if (length(parm) > 1) {
            for (l in parm) coda::densplot(object$mcmc[[l]], ...)
        } else {
            coda::densplot(object$mcmc[[parm]], ...)
        }
    }
    invisible()
}

cumuplot <- function (object, ...) UseMethod("cumuplot")

cumuplot.jm <- function (object,
                         parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                                  "tau_bs_gammas", "gammas", "alphas"), ...) {
    parm <- match.arg(parm)
    if (parm == "all") {
        nams_parms <- c("betas", "sigmas", "D", "bs_gammas", "tau_bs_gammas",
                        "gammas", "alphas")
        nams_mcmc <- names(object$mcmc)
        ind <- unlist(sapply(paste0("^", nams_parms), grep, nams_mcmc))
        nams_mcmc <- nams_mcmc[ind]
        for (i in seq_along(nams_mcmc)) {
            parms_i <- nams_mcmc[[i]]
            x <- object$mcmc[[parms_i]]
            if (!is.null(x)) coda::cumuplot(x, ...)
        }
    } else {
        parm <- grep(paste0("^", parm), names(object$mcmc))
        if (length(parm) > 1) {
            for (l in parm) coda::cumuplot(object$mcmc[[l]], ...)
        } else {
            coda::cumuplot(object$mcmc[[parm]], ...)
        }
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
                          `97.5%` = object$statistics$CI_upp[[name]],
                          P = object$statistics$P[[name]],
                          row.names = names(object$statistics$P[[name]]),
                          check.names = FALSE)
        Rhat <- object$statistics$Rhat[[name]][, 1L]
        if (!is.null(Rhat))
            out$Rhat <- Rhat
        out
    }
    fam_names <- sapply(families, "[[", "family")
    has_sigma_fam <- c("gaussian", "Student-t", "beta", "Gamma",
                   "negative binomial", "beta binomial", "censored normal")
    has_sigmas <- object$model_data$has_sigmas
    has_sigmas[has_sigmas > 0] <- which(has_sigmas > 0)
    tab_sigmas <- tab_f("sigmas")
    for (i in seq_len(n_outcomes)) {
        nam_outcome <- paste0("Outcome", i)
        out[[nam_outcome]] <- tab_f(paste0("betas", i))
        if (fam_names[i] %in% has_sigma_fam) {
            k <- nrow(out[[nam_outcome]])
            out[[nam_outcome]] <-
                rbind(out[[nam_outcome]],
                      tab_sigmas[paste0("sigmas_", has_sigmas[i]), ])
            row.names(out[[nam_outcome]])[k + 1] <- "sigma"
        }
    }
    out$Survival <- do.call(rbind, list(tab_f("gammas"), tab_f("alphas")))
    out$fit_stats <- object$fit_stats
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
    if (!is.null(x$fit_stats$conditional$DIC)) {
        model.sum <-
            data.frame(DIC = c(x$fit_stats$marginal$DIC, x$fit_stats$conditional$DIC),
                       WAIC = c(x$fit_stats$marginal$WAIC, x$fit_stats$conditional$WAIC),
                       LPML = c(x$fit_stats$marginal$LPML, x$fit_stats$conditional$LPML),
                       row.names = c("marginal", "conditional"))
        cat("\n")
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
        "\nthinning:", x$control$n_thin,
        "\ntime:", if (tt < 1) {round(tt * 60)} else if (tt > 60)
            {round(tt/60, 1)} else {round(tt, 1)},
        if (tt < 1) {"sec"} else if (tt > 60) {"hours"} else {"min"})
    cat("\n")
    invisible(x)
}

print.jm <- function (x, digits = max(4, getOption("digits") - 4), ...) {
    xx <- summary(x)
    cat("\nCall:\n", printCall(xx$call), "\n", sep = "")
    cat("\nRandom-effects covariance matrix:\n")
    D <- xx$D
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

coef.jm <- function (object, ...) {
    gammas <- object$statistics$Mean[["gammas"]]
    if (is.null(gammas)) object$statistics$Mean[["alphas"]] else
        list("gammas" = gammas,
             "association" = object$statistics$Mean[["alphas"]])
}

fixef.jm <- function(object, outcome = Inf, ...) {
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

ranef.jm <- function(object, outcome = Inf, post_vars = FALSE, ...) {
    if (!is.numeric(outcome) || outcome < 0) {
        stop("'outcome' should be a positive integer.")
    }
    outcome <- round(outcome)
    if (outcome > length(object$model_data$y)) {
        out <- object$statistics$Mean$b
        if (post_vars)
            attr(out, "post_vars") <- object$statistics$post_vars
    } else {
        ind <- object$model_data$ind_RE[[outcome]]
        out <- object$statistics$Mean$b[, ind, drop = FALSE]
        if (post_vars)
            attr(out, "post_vars") <-
            object$statistics$post_vars[ind, ind, , drop = FALSE]
    }
    out
}

terms.jm <- function (x, process = c("longitudinal", "event"),
                      type = c("fixed", "random"), ...) {
    process <- match.arg(process)
    type <- match.arg(type)
    combo <- paste(process, type, sep = "_")
    switch(combo,
           "longitudinal_fixed" = x$model_info$terms$terms_FE,
           "longitudinal_random" = x$model_info$terms$terms_RE,
           "event_fixed" = , "event_random" = x$model_info$terms$terms_Surv)
}

model.frame.jm <- function (formula, process = c("longitudinal", "event"),
                            type = c("fixed", "random"), ...) {
    process <- match.arg(process)
    type <- match.arg(type)
    combo <- paste(process, type, sep = "_")
    switch(combo,
           "longitudinal_fixed" = formula$model_info$frames$mf_FE,
           "longitudinal_random" = formula$model_info$frames$mf_RE,
           "event_fixed" = , "event_random" = formula$model_info$frames$mf_Surv)
}

model.matrix.jm <- function (object, ...) {
    tr <- terms(object)
    mf <- model.frame(object)
    if (is.data.frame(mf)) {
        model.matrix(tr, mf)
    } else {
        mapply(model.matrix.default, object = tr, data = mf, SIMPLIFY = FALSE)
    }
}

family.jm <- function (object, ...) {
    object$model_info$families
}

ggtraceplot <- function (object, ...) UseMethod("ggtraceplot")

ggtraceplot.jm <- function(object,
                        parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                                 "tau_bs_gammas", "gammas", "alphas"),
                        size = 1, alpha = 0.8,
                        theme = c('standard', 'catalog', 'metro',
                                      'pastel', 'beach', 'moonlight', 'goo',
                                      'sunset'), grid = FALSE,
                        gridrows = 3, gridcols = 1,
                        ...) {
    chain <- iteration <- NULL
    parm <- match.arg(parm)
    coltheme <- match.arg(theme)
    ggdata <- ggprepare(object, parm)
    n_parms <- length(unique(ggdata$parm))
    n_chains <- object$control$n_chains
    if (grid) {
        gplots <- list(NULL)
        for (i in seq_len(n_parms)) {
            gplots[[i]] <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
                geom_line(aes(x = iteration, y = value, color = chain),
                          size = size, alpha = alpha) +
                ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = ggcolthemes[[coltheme]]) +
                guides(color = guide_legend(override.aes = list(alpha = 1)))
        }
        marrangeGrob(grobs = gplots, nrow = gridrows, ncol = gridcols)
    } else {
        for (i in seq_len(n_parms)) {
            g <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
                geom_line(aes(x = iteration, y = value, color = chain),
                          size = size, alpha = alpha) +
                ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = ggcolthemes[[coltheme]]) +
                guides(color = guide_legend(override.aes = list(alpha = 1)))
            print(g)
        }
    }
}

ggdensityplot <- function (object, ...) UseMethod("ggdensityplot")

ggdensityplot.jm <- function(object,
                      parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                               "tau_bs_gammas", "gammas", "alphas"),
                      size = 1, alpha = 0.6,
                      theme = c('standard', 'catalog', 'metro',
                                    'pastel', 'beach', 'moonlight', 'goo',
                                    'sunset'), grid = FALSE,
                      gridrows = 3, gridcols = 1,
                      ...) {
    chain <- NULL
    parm <- match.arg(parm)
    coltheme <- match.arg(theme)
    ggdata <- ggprepare(object, parm)
    n_parms <- length(unique(ggdata$parm))
    n_chains <- object$control$n_chains
    if (grid) {
        gplots <- list(NULL)
        for (i in seq_len(n_parms)) {
            gplots[[i]] <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
                geom_density(aes(x = value, color = chain, fill = chain),
                             size = size, alpha = alpha) +
                ggtitle(paste('Density plot of ', unique(ggdata$parm)[i])) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = ggcolthemes[[coltheme]]) +
                scale_fill_manual(values = ggcolthemes[[coltheme]]) +
                guides(color = guide_legend(override.aes = list(alpha = 1)))
        }
        marrangeGrob(grobs = gplots, nrow = gridrows, ncol = gridcols)
    } else {
        for (i in seq_len(n_parms)) {
            g <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
                geom_density(aes(x = value, color = chain, fill = chain),
                             size = size, alpha = alpha) +
                ggtitle(paste('Density plot of ', unique(ggdata$parm)[i])) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = ggcolthemes[[coltheme]]) +
                scale_fill_manual(values = ggcolthemes[[coltheme]]) +
                guides(color = guide_legend(override.aes = list(alpha = 1)))
            print(g)
        }
    }
}

effectPlotData.jm <- function (object, newdata, level = 0.95, ...) {
    termsX <- object$model_info$terms$terms_FE_noResp
    xlevels <- mapply2(.getXlevels, termsX, object$model_info$frames$mf_FE)
    mfX <- mapply2(model.frame.default, formula = termsX, xlev = xlevels,
                   MoreArgs = list(data = newdata))
    X <- mapply2(model.matrix.default, object = termsX, data = mfX)
    ind_betas <- grep("^betas", names(object$mcmc))
    betas <- object$mcmc[ind_betas]
    betas <- lapply(betas, function (b) do.call("rbind", b))
    Xbetas <- mapply2(tcrossprod, X, betas)
    pred <- lapply(Xbetas, rowMeans)
    names(pred) <- paste0("pred", seq_along(pred))
    Qs <- lapply(Xbetas, rowQuantiles,
                 probs = c((1 - level) / 2, (1 + level) / 2))
    for (i in seq_along(Qs)) {
        colnames(Qs[[i]]) <- paste0(c("low", "upp"), i)
    }
    cbind(newdata, do.call("cbind", pred), do.call("cbind", Qs))
}

compare_jm <- function (..., type = c("marginal", "conditional"),
                        order = c("WAIC", "DIC", "LPML", "none")) {
    model_names <- sapply(substitute(list(...)), deparse)[-1L]
    models <- list(...)
    if (!all(sapply(models, inherits, "jm"))) {
        stop("compare_jm() works with jm objects.")
    }
    if (length(models) == 1L) {
        stop("compare_jm() is supposed to compare two or more joint models.")
    }
    #respVars <- lapply(models, function (m) m$model_info$var_names$respVars)
    #check_names <- sapply(respVars[-1],
    #                      function (nams, nams_1) all(nams %in% nams_1),
    #                      nams_1 = respVars[[1]])
    #if (!all(check_names)) {
    #    stop("it seems that some joint have different longitudinal outcomes.")
    #}
    type <- match.arg(type)
    order <- match.arg(order)
    extract_criteria <- function (m, type) {
        if (type == "marginal") {
            data.frame(DIC = m$fit_stats$marginal$DIC,
                       WAIC = m$fit_stats$marginal$WAIC,
                       LPML = m$fit_stats$marginal$LPML, check.names = FALSE)
        } else {
            data.frame(DIC = m$fit_stats$conditional$DIC,
                       WAIC = m$fit_stats$conditional$WAIC,
                       LPML = m$fit_stats$conditional$LPML, check.names = FALSE)
        }
    }
    out <- do.call("rbind", lapply(models, extract_criteria, type = type))
    out$model <- model_names
    out <- out[c("model", "DIC", "WAIC", "LPML")]
    names(out) <- c(" ", "DIC", "WAIC", "LPML")
    if (order != "none") {
        out <- if (order == "LPML") out[order(out[[order]], decreasing = TRUE), ]
            else out[order(out[[order]]), ]
    }
    out <- list(table = out, type = type)
    class(out) <- "compare_jm"
    out
}

print.compare_jm <- function (x, ...) {
    cat("\n")
    print.data.frame(x$table, row.names = FALSE)
    cat("\nThe criteria are calculated based on the", x$type, "log-likelihood.")
}

crLong <- function (data, statusVar, censLevel, nameStrata = "strata",
          nameStatus = "status2") {
    n <- nrow(data)
    status <- data[[statusVar]]
    unqLevs <- unique(status)
    unqLevs <- unqLevs[unqLevs != censLevel]
    ncr <- length(unqLevs)
    dataOut <- data[rep(seq_len(n), each = ncr), , drop = FALSE]
    dataOut[[nameStrata]] <- rep(unqLevs, n)
    dataOut[[nameStatus]] <- as.numeric(dataOut[[statusVar]] ==
                                            dataOut[[nameStrata]])
    dataOut[[nameStrata]] <- factor(dataOut[[nameStrata]])
    dataOut
}

predict.jm <- function (object, newdata = NULL, newdata2 = NULL,
                        times = NULL,
                        process = c("longitudinal", "event"),
                        type_pred = c("response", "link"),
                        type = c("subject_specific", "mean_subject"),
                        level = 0.95, return_newdata = FALSE,
                        n_samples = 200L, n_mcmc = 55L, cores = NULL,
                        seed = 123L, ...) {
    process <- match.arg(process)
    type_pred <- match.arg(type_pred)
    type <- match.arg(type)
    id_var <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    Time_var <- object$model_info$var_names$Time_var
    event_var <- object$model_info$var_names$event_var
    respVars <- unlist(object$model_info$var_names$respVars)
    if (object$model_info$CR_MS && is.data.frame(newdata)) {
        stop("for competing risks and multi-state models, argument 'newdata' ",
             "must be a list of two data.frames, one for the longitudinal ",
             "outcomes and one for the event process, the latter under the ",
             "correct format.\n")
    }
    if (!is.data.frame(newdata)) {
        if (!is.list(newdata) || length(newdata) != 2
            || !names(newdata) %in% c("newdataL", "newdataE")) {
            stop("'newdata' must be a list with two data.frame elements ",
                 "named 'newdataL' and 'newdataE'.\n")
        }
        for (i in seq_along(respVars)) {
            v <- respVars[i]
            if (is.null(newdata$newdataE[[v]])) {
                newdata$newdataE[[v]] <- rep(0.1, nrow(newdata$newdataE))
            }
        }
        termsL <- object$model_info$terms$terms_FE_noResp
        all_vars <- unlist(lapply(termsL, all.vars), use.names = FALSE)
        all_vars <- all_vars[!all_vars %in% time_var]
        missing_vars <- all_vars[!all_vars %in% names(newdata$newdataE)]
        if (length(missing_vars)) {
            stop("the data.frame 'newdata$newdataE' should contain the ",
                 "variable(s): ", paste(missing_vars, collapse = ", "), ".\n")
        }
        missing_vars <- all_vars[!all_vars %in% names(newdata$newdataL)]
        if (length(missing_vars)) {
            stop("the data.frame 'newdata$newdataL' should contain the ",
                 "variable(s): ", paste(missing_vars, collapse = ", "), ".\n")
        }
    }
    if (is.data.frame(newdata)) {
        if (is.null(newdata[[event_var]])) newdata[[event_var]] <- 0
        if (is.null(newdata[[Time_var]])) {
            last_time <- function (x) max(x, na.rm = TRUE)
            f <- factor(newdata[[id_var]], unique(newdata[[id_var]]))
            newdata[[Time_var]] <- ave(newdata[[time_var]], f, FUN = last_time)
        }
        termsL <- object$model_info$terms$terms_FE_noResp
        all_vars <- unlist(lapply(termsL, all.vars), use.names = FALSE)
        all_vars <- all_vars[!all_vars %in% time_var]
        all_vars <- c(all_vars, all.vars(object$model_info$terms$terms_Surv_noResp))
        missing_vars <- all_vars[!all_vars %in% names(newdata)]
        if (length(missing_vars)) {
            stop("the data.frame 'newdata' should contain the ",
                 "variable(s): ", paste(missing_vars, collapse = ", "), ".\n")
        }
    }
    if (is.null(cores)) {
        n <- if (!is.data.frame(newdata)) length(unique(newdata$newdataL[[id_var]]))
        else length(unique(newdata[[id_var]]))
        cores <- if (n > 20) 4L else 1L
    }
    components_newdata <- get_components_newdata(object, newdata, n_samples,
                                                 n_mcmc, cores, seed)
    if (process == "longitudinal") {
        predict_Long(object, components_newdata, newdata, newdata2, times, type,
                     type_pred, level, return_newdata)
    } else {
        predict_Event(object, components_newdata, newdata, newdata2, times,
                      level, return_newdata)
    }
}

plot.predict_jm <- function (x, x2 = NULL, subject = 1, outcomes = 1,
                             fun_long = NULL, fun_event = NULL,
                             CI_long = TRUE, CI_event = TRUE,
                             xlab = "Follow-up Time", ylab_long = NULL,
                             ylab_event = "Cumulative Risk", main = "",
                             lwd_long = 2, lwd_event = 2,
                             ylim_long_outcome_range = TRUE,
                             col_line_long = "#0000FF",
                             col_line_event = c("#FF0000", "#03BF3D", "#8000FF"),
                             pch_points = 16, col_points = "blue", cex_points = 1,
                             fill_CI_long = "#0000FF4D",
                             fill_CI_event = c("#FF00004D", "#03BF3D4D", "#8000FF4D"),
                             cex_xlab = 1, cex_ylab_long = 1, cex_ylab_event = 1,
                             cex_main = 1, cex_axis = 1, col_axis = "black",
                             pos_ylab_long = c(0.1, 2, 0.08), bg = "white",
                             ...) {
    process_x <- attr(x, "process")
    pred_Long <- if (process_x == "longitudinal") x
    pred_Event <- if (process_x == "event") x
    if (!is.null(x2)) {
        process_x2 <- attr(x2, "process")
        if (process_x2 == "longitudinal" && is.null(pred_Long)) pred_Long <- x2
        if (process_x2 == "event" && is.null(pred_Event)) pred_Event <- x2
    }
    id_var <- attr(x, "id_var")
    time_var <- attr(x, "time_var")
    resp_vars <- attr(x, "resp_vars")
    ranges <- attr(x, "ranges")
    last_times <- attr(x, "last_times")
    y <- attr(x, "y")
    times_y <- attr(x, "times_y")
    id <- attr(x, "id")
    if (!is.null(pred_Long)) {
        test1 <- is.data.frame(pred_Long)
        test2 <- is.list(pred_Long) && length(pred_Long) == 2L && is.data.frame(pred_Long[[1]])
        if (!test1 && !test2) {
            stop("you must use predict.jm(..., return_newdata = TRUE)")
        }
        if (test2) {
            pred_Long <- rbind(pred_Long[[1L]], pred_Long[[2L]])
        }
    }
    if (!is.null(pred_Event) && !is.data.frame(pred_Event)) {
        stop("you must use predict.jm(..., return_newdata = TRUE)")
    }
    unq_id <- if (!is.null(pred_Long)) unique(pred_Long[[id_var]])
    if (!is.null(pred_Event)) unq_id <- unique(c(pred_Event[[id_var]], unq_id))
    if (length(subject) > 1L) {
        stop("'subject' must be of length 1.")
    }
    if (!subject %in% unq_id && subject > length(unq_id)) {
        stop("not valid input for 'subject'.")
    }
    subj <- if (subject %in% unq_id) subject else unq_id[subject]
    subj_ind <- match(subj, unq_id)
    if (!is.null(pred_Long)) {
        pred_Long <- pred_Long[pred_Long[[id_var]] == subj, ]
        if (!is.null(pred_Event)) {
            pred_Long <- pred_Long[pred_Long[[time_var]] <= last_times[subj_ind], ]
        }
        if (!nrow(pred_Long)) {
            stop("no available measurements before the last time.")
        }
        pos_outcomes <- grep("pred_", names(pred_Long), fixed = TRUE)
        n_outcomes <- length(outcomes)
        if (n_outcomes > length(pos_outcomes)) {
            stop("the length of 'outcomes' is greater than the number of ",
                 "outcomes in the dataset.")
        }
        if (any(outcomes > length(pos_outcomes))) {
            stop("not valid entries in 'outcome'.")
        }
        if (!is.null(pred_Event) && n_outcomes > 3) {
            warning("when 'pred_Event' is not null max three outcomes are allowed in the plot.")
            n_outcomes <- 3
            outcomes <- rep_len(outcomes, length.out = 3L)
        }
        if (is.null(fun_long)) {
            fun_long <- rep(list(function (x) x), n_outcomes)
        } else {
            if (is.function(fun_long)) fun_long <- rep(list(fun_long), n_outcomes)
            if (is.list(fun_long) && (length(fun_long) != n_outcomes ||
                                      !all(sapply(fun_long, is.function)))) {
                stop("'fun_long' needs to be a function or a list of functions.")
            }
        }
        col_line_long <- rep(col_line_long, length.out = n_outcomes)
        pch_points <- rep(pch_points, length.out = n_outcomes)
        col_points <- rep(col_points, length.out = n_outcomes)
        cex_points <- rep(cex_points, length.out = n_outcomes)
        fill_CI_long <- rep(fill_CI_long, length.out = n_outcomes)
    }
    if (!is.null(pred_Event)) {
        pred_Event <- pred_Event[pred_Event[[id_var]] == subj, ]
        if (is.null(fun_event) || !is.function(fun_event)) {
            fun_event <- function (x) x
        }
    }
    if (is.null(ylab_long)) {
        ylab_long <- resp_vars
    }
    xlim <- NULL
    if (!is.null(pred_Long)) xlim <- range(xlim, pred_Long[[time_var]])
    if (!is.null(pred_Event)) xlim <- range(xlim, pred_Event[[time_var]])
    plot_long_i <- function (outcome, add_xlab = FALSE, box = TRUE,
                             cex_axis = cex_axis) {
        ind <- pos_outcomes[outcome]
        outcome_i <- match(outcome, outcomes)
        f <- fun_long[[outcome_i]]
        preds <- f(pred_Long[[ind]])
        low <- f(pred_Long[[ind + 1]])
        upp <- f(pred_Long[[ind + 2]])
        times <- pred_Long[[time_var]]
        ry <- range(preds, low, upp)
        rx <- range(times)
        y_lim <- if (ylim_long_outcome_range) {
            range(f(ranges[[outcome]]), ry)
        } else {
            ry
        }
        plot(rx, ry, type = "n", xaxt = "n", bty = if (box) "o" else "n",
             xlab = if (add_xlab) xlab  else "", xlim = xlim, col.axis = col_axis,
             ylim = y_lim, ylab = ylab_long[outcome],
             cex.lab = cex_ylab_long, cex.axis = cex_axis, col.lab = col_axis,
             col.axis = col_axis)
        if (!add_xlab) {
            axis(1, c(-5, last_times[subj_ind]), labels = c("", ""), tcl = 0,
                 cex.axis = cex_axis, col = col_axis, col.axis = col_axis,
                 col.ticks = col_axis)
        }
        if (CI_long) {
            polygon(c(times, rev(times)), c(low, rev(upp)), border = NA,
                    col = fill_CI_long[outcome_i])
        }
        y_i <- f(c(y[[outcome]]))
        times_y_i <- times_y[[outcome]]
        id_i <- id[[outcome]]
        points(times_y_i[id_i == subj_ind], y_i[id_i == subj_ind],
               pch = pch_points[outcome_i], cex = cex_points[outcome_i],
               col = col_points[outcome_i])
        lines(times, preds, lwd = lwd_long, col = col_line_long[outcome_i])
        abline(v = last_times[subj_ind] + 0.01, lty = 3, col = col_axis)
    }
    plot_event <- function (box = FALSE, axis_side = 4, cex_axis = cex_axis) {
        ind <- grep("pred_", names(pred_Event), fixed = TRUE)
        preds <- fun_event(pred_Event[[ind]])
        low <- fun_event(pred_Event[[ind + 1]])
        upp <- fun_event(pred_Event[[ind + 2]])
        strata <- pred_Event[["_strata"]]
        if (is.null(strata)) strata <- rep(1, length(preds))
        unq_strata <- unique(strata)
        col_line_event <- rep(col_line_event, length.out = length(unq_strata))
        fill_CI_event <- rep(fill_CI_event, length.out = length(unq_strata))
        times <- pred_Event[[time_var]]
        ry <- sort(fun_event(c(0, 1)))
        rx <- range(times)
        plot(rx, ry, type = "n", xlab = "", ylab = "", xlim = xlim,
             axes = FALSE, col.axis = col_axis, col.lab = col_axis, ylim = ry)
        if (box) box(col = col_axis)
        axis(axis_side, cex.axis = cex_axis, col = col_axis,
             col.ticks = col_axis, col.axis = col_axis)
        for (i in seq_along(unq_strata)) {
            ind_str <- strata == unq_strata[i]
            if (CI_event) {
                polygon(c(times[ind_str], rev(times[ind_str])),
                        c(low[ind_str], rev(upp[ind_str])), border = NA,
                        col = fill_CI_event[i])
            }
            lines(times[ind_str], preds[ind_str], lwd = lwd_event,
                  col = col_line_event[i])
        }
    }
    if (is.null(pred_Event)) {
        for (i in seq_along(outcomes)) {
            plot_long_i(outcomes[i], TRUE, cex_axis = cex_axis)
            title(main = main, cex = cex_main)
            axis(1, cex.axis = cex_axis, col = col_axis,
                 col.ticks = col_axis, col.axis = col_axis)
        }
    }
    if (is.null(pred_Long)) {
        plot_event(box = TRUE, 2, cex_axis = cex_axis)
        title(xlab = xlab, cex = cex_xlab)
        title(ylab = ylab_event, cex = cex_ylab_event)
        title(main = main, cex = cex_main)
        abline(v = last_times[subj_ind] + 0.01, lty = 3)
        axis(1, cex.axis = cex_axis, col = col_axis, col.ticks = col_axis,
             col.axis = col_axis)
    }
    if (!is.null(pred_Long) && !is.null(pred_Event)) {
        if (n_outcomes == 1) {
            # n_outcomes == 1
            op <- par(mar = c(4,4,3,4), mgp = c(2, 0.4, 0), tcl = -0.3, bg = bg)
            plot_long_i(outcomes[1L], cex_axis = cex_axis)
            axis(1, cex.axis = cex_axis, col.ticks = col_axis, col = col_axis,
                 col.axis = col_axis)
            title(xlab = xlab, cex = cex_xlab, col = col_axis)
            par(new = TRUE)
            plot_event(cex_axis = cex_axis)
            mtext(ylab_event, 4, 1.5, cex = cex_ylab_event, col = col_axis)
            par(op)
            mtext(main, 3, 1.5, cex = cex_main, col = col_axis)
        } else if (n_outcomes == 2) {
            # n_outcomes == 2
            op <- par(mfrow = c(2, 1), oma = c(4,4,3,4), mar = c(0, 0, 0, 0),
                      mgp = c(2, 0.4, 0), tcl = -0.3, bg = bg)
            pp <- par("usr")[1] + pos_ylab_long * diff(par("usr")[1:2])
            plot_long_i(outcomes[1L], box = FALSE, cex_axis = cex_axis)
            mtext(ylab_long[outcomes[1L]], 2, 1.5, at = pp[1],
                  cex = cex_ylab_long * 0.66, col = col_axis)
            plot_long_i(outcomes[2L], box = FALSE, cex_axis = cex_axis)
            mtext(ylab_long[outcomes[2L]], 2, 1.5, at = pp[2],
                  cex = cex_ylab_long * 0.66, col = col_axis)
            axis(1, cex.axis = cex_axis, col.ticks = col_axis, col = col_axis,
                 col.axis = col_axis)
            mtext(xlab, side = 1, line = 1.5, outer = TRUE,
                  cex = cex_xlab, col = col_axis)
            par(op)
            op <- par(new = TRUE, oma = c(4,4,3,4), mar = c(0, 0, 0, 0),
                      mgp = c(2, 0.4, 0), tcl = -0.3, cex = 0.9)
            plot_event(box = TRUE, cex_axis = 0.66 * cex_axis)
            mtext(ylab_event, 4, 1.5, cex = cex_ylab_event, col = col_axis)
            par(op)
            mtext(main, 3, 1.5, cex = cex_main, col = col_axis)
        } else {
            # n_outcomes == 3
            op <- par(mfrow = c(3, 1), oma = c(4,4,3,4), mar = c(0, 0, 0, 0),
                      mgp = c(2, 0.4, 0), tcl = -0.3, bg = bg)
            pp <- par("usr")[1] + pos_ylab_long * diff(par("usr")[1:2])
            plot_long_i(outcomes[1L], box = FALSE, cex_axis = cex_axis)
            mtext(ylab_long[outcomes[1L]], 2, 1.5, at = pp[1],
                  cex = cex_ylab_long * 0.66, col = col_axis)
            plot_long_i(outcomes[2L], box = FALSE, cex_axis = cex_axis)
            mtext(ylab_long[outcomes[2L]], 2, 1.5, at = pp[2],
                  cex = cex_ylab_long * 0.66, col = col_axis)
            plot_long_i(outcomes[3L], box = FALSE, cex_axis = cex_axis)
            mtext(ylab_long[outcomes[3L]], 2, 1.5, at = pp[3],
                  cex = cex_ylab_long * 0.66, col = col_axis)
            axis(1, cex.axis = cex_axis, col = col_axis, col.ticks = col_axis,
                 col.axis = col_axis)
            mtext(xlab, side = 1, line = 1.5, outer = TRUE, cex = cex_xlab,
                  col = col_axis)
            box("inner", col = col_axis)
            par(op)
            op <- par(new = TRUE, oma = 0.6525 * c(4,4,3,4), mar = c(0, 0, 0, 0),
                      mgp = c(2, 0.4, 0), tcl = -0.3, cex = 0.66)
            plot_event(cex_axis = cex_axis)
            mtext(ylab_event, 4, 1.5, cex = cex_ylab_event, col = col_axis)
            par(op)
            mtext(main, 3, 1.5, cex = cex_main, col = col_axis)
        }
    }
    invisible()
}

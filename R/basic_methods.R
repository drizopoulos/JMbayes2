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
        ind <- unlist(sapply(paste0("^", nams_parms), grep, nams_mcmc))
        nams_mcmc <- nams_mcmc[ind]
        out <- vector("list", length(nams_mcmc))
        names(out) <- nams_mcmc
        for (i in seq_along(out)) {
            parms_i <- nams_mcmc[[i]]
            x <- object$mcmc[[parms_i]]
            if (!is.null(x)) out[[i]] <- coda::gelman.diag(x)
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

densplot <- function (object, ...) UseMethod("densityplot")

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
        ind <- object$model_data$ind_RE[outcome]
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

model.matrix.jm <- function (object, process = c("longitudinal", "event"),
                             type = c("fixed", "random"), ...) {
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

# ggplot mcmc diagnostics need ggplot2 and gridExtra
ggtraceplot <- function (object, ...) UseMethod("ggtraceplot")

# traceplot with ggplot
ggtraceplot.jm <- function(object,
                        parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                                 "tau_bs_gammas", "gammas", "alphas"),
                        size = 1, alpha = 0.8,
                        chaincols = c('standard', 'catalog', 'metro',
                                      'pastel', 'beach', 'moonlight', 'goo',
                                      'sunset'),
                        gridrow = 3, gridcol = 1, grid = FALSE,
                        ...) {
    chain <- iteration <- NULL
    parm <- match.arg(parm)
    coltheme <- match.arg(chaincols)
    ggdata <- ggprepare(object, parm)
    n_parms <- length(unique(ggdata$parm))
    n_chains <- object$control$n_chains
    if (grid) {
        gplots <- list(NULL)
        for (i in seq_len(n_parms)) {
            gplots[[i]] <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
                geom_line(aes(x = iteration, y = value, color = chain), size = size, alpha = alpha) +
                ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = ggcolthemes[[coltheme]]) +
                guides(color = guide_legend(override.aes = list(alpha = 1)))
        }
        marrangeGrob(grobs = gplots, nrow = gridrow, ncol = gridcol)
    } else {
        for (i in seq_len(n_parms)) {
            g <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
                geom_line(aes(x = iteration, y = value, color = chain), size = size, alpha = alpha) +
                ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = ggcolthemes[[coltheme]]) +
                guides(color = guide_legend(override.aes = list(alpha = 1)))
            print(g)
        }
    }
}

ggdensityplot <- function (object, ...) UseMethod("ggdensityplot")

# density plot with ggplot
ggdensityplot.jm <- function(object,
                      parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                               "tau_bs_gammas", "gammas", "alphas"),
                      size = 1, alpha = 0.6,
                      chaincols = c('standard', 'catalog', 'metro',
                                    'pastel', 'beach', 'moonlight', 'goo',
                                    'sunset'),
                      gridrow = 3, gridcol = 1, grid = FALSE,
                      ...) {
    chain <- NULL
    parm <- match.arg(parm)
    coltheme <- match.arg(chaincols)
    ggdata <- ggprepare(object, parm)
    n_parms <- length(unique(ggdata$parm))
    n_chains <- object$control$n_chains
    if (grid) {
        gplots <- list(NULL)
        for (i in seq_len(n_parms)) {
            gplots[[i]] <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
                geom_density(aes(x = value, color = chain, fill = chain), size = size, alpha = alpha) +
                ggtitle(paste('Density plot of ', unique(ggdata$parm)[i])) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = ggcolthemes[[coltheme]]) +
                scale_fill_manual(values = ggcolthemes[[coltheme]]) +
                guides(color = guide_legend(override.aes = list(alpha = 1)))
        }
        marrangeGrob(grobs = gplots, nrow = gridrow, ncol = gridcol)
    } else {
        for (i in seq_len(n_parms)) {
            g <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
                geom_density(aes(x = value, color = chain, fill = chain), size = size, alpha = alpha) +
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
    respVars <- lapply(models, function (m) m$model_info$var_names$respVars)
    check_names <- sapply(respVars[-1],
                          function (nams, nams_1) all(nams %in% nams_1),
                          nams_1 = respVars[[1]])
    if (!all(check_names)) {
        stop("it seems that some joint have different longitudinal outcomes.")
    }
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




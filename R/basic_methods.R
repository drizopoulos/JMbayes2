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
                call = object$call, recurrent = object$model_info$recurrent,
                any_terminal = length(object$model_data$which_term_h) > 0)
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
    if(out$recurrent & out$any_terminal) {
      out$Survival <- do.call(rbind, list(tab_f("gammas"), tab_f("alphas"),
                                          tab_f("alphaF")))
    } else {
      out$Survival <- do.call(rbind, list(tab_f("gammas"), tab_f("alphas")))
    }
    out$sigmaF <- tab_f("sigmaF")[c(1, 3, 4)]
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
    if(x$recurrent) {
      cat("\nFrailty standard deviation:\n")
      print(round(x[["sigmaF"]], digits))
    }
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
    if(xx$recurrent) {
      cat("\nFrailty standard deviation:\n")
      print(round(xx[["sigmaF"]], digits))
    }
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
                                     'sunset', 'custom'), grid = FALSE,
                           gridrows = 3, gridcols = 1, custom_theme = NULL,
                           ...) {
  chain <- iteration <- NULL
  parm <- match.arg(parm)
  coltheme <- match.arg(theme)
  ggdata <- ggprepare(object, parm)
  n_parms <- length(unique(ggdata$parm))
  n_chains <- object$control$n_chains
  if(!is.null(custom_theme)) {
    if (length(custom_theme) != n_chains)
      stop('User specified custom color themes should be a named character vector with one color specified for each chain')
    coltheme <- 'custom'
    ggcolthemes[[coltheme]] <- custom_theme
  }
  if (grid) {
    gplots <- list(NULL)
    for (i in seq_len(n_parms)) {
      if (n_chains == 3 | !is.null(custom_theme)) {
        gplots[[i]] <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
          geom_line(aes(x = iteration, y = value, color = chain),
                    size = size, alpha = alpha) +
          ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
          theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
          scale_color_manual(values = ggcolthemes[[coltheme]]) +
          guides(color = guide_legend(override.aes = list(alpha = 1)))
      } else {
        gplots[[i]] <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
          geom_line(aes(x = iteration, y = value, color = chain),
                    size = size, alpha = alpha) +
          ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
          theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
          guides(color = guide_legend(override.aes = list(alpha = 1)))
      }
    }
    marrangeGrob(grobs = gplots, nrow = gridrows, ncol = gridcols)
  } else {
    for (i in seq_len(n_parms)) {
      if (n_chains == 3 | !is.null(custom_theme)) {
        g <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
          geom_line(aes(x = iteration, y = value, color = chain),
                    size = size, alpha = alpha) +
          ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
          theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
          scale_color_manual(values = ggcolthemes[[coltheme]]) +
          guides(color = guide_legend(override.aes = list(alpha = 1)))
        print(g)
      } else {
        g <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
          geom_line(aes(x = iteration, y = value, color = chain),
                    size = size, alpha = alpha) +
          ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
          theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
          guides(color = guide_legend(override.aes = list(alpha = 1)))
        print(g)
      }
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
                                'sunset', 'custom'), grid = FALSE,
                      gridrows = 3, gridcols = 1, custom_theme = NULL,
                      ...) {
    chain <- NULL
    parm <- match.arg(parm)
    coltheme <- match.arg(theme)
    ggdata <- ggprepare(object, parm)
    n_parms <- length(unique(ggdata$parm))
    n_chains <- object$control$n_chains
    if(!is.null(custom_theme)) {
      if (length(custom_theme) != n_chains)
        stop('User specified custom color themes should be a named character vector with one color specified for each chain')
      coltheme <- 'custom'
      ggcolthemes[[coltheme]] <- custom_theme
    }
    if (grid) {
        gplots <- list(NULL)
        for (i in seq_len(n_parms)) {
          if (n_chains == 3 | !is.null(custom_theme)) {
            gplots[[i]] <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
              geom_density(aes(x = value, color = chain, fill = chain),
                           size = size, alpha = alpha) +
              ggtitle(paste('Density plot of ', unique(ggdata$parm)[i])) +
              theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
              scale_color_manual(values = ggcolthemes[[coltheme]]) +
              scale_fill_manual(values = ggcolthemes[[coltheme]]) +
              guides(color = guide_legend(override.aes = list(alpha = 1)))
          } else {
            gplots[[i]] <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
              geom_density(aes(x = value, color = chain, fill = chain),
                           size = size, alpha = alpha) +
              ggtitle(paste('Density plot of ', unique(ggdata$parm)[i])) +
              theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
              guides(color = guide_legend(override.aes = list(alpha = 1)))
          }
        }
        marrangeGrob(grobs = gplots, nrow = gridrows, ncol = gridcols)
    } else {
        for (i in seq_len(n_parms)) {
          if (n_chains == 3 | !is.null(custom_theme)) {
            g <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
              geom_density(aes(x = value, color = chain, fill = chain),
                           size = size, alpha = alpha) +
              ggtitle(paste('Density plot of ', unique(ggdata$parm)[i])) +
              theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
              scale_color_manual(values = ggcolthemes[[coltheme]]) +
              scale_fill_manual(values = ggcolthemes[[coltheme]]) +
              guides(color = guide_legend(override.aes = list(alpha = 1)))
            print(g)
          } else {
            g <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
              geom_density(aes(x = value, color = chain, fill = chain),
                           size = size, alpha = alpha) +
              ggtitle(paste('Density plot of ', unique(ggdata$parm)[i])) +
              theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
              guides(color = guide_legend(override.aes = list(alpha = 1)))
            print(g)
          }
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

crisk_setup <- function (data, statusVar, censLevel, nameStrata = "strata",
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

predict.jm <- function (object, newdata = NULL, newdata2 = NULL, times = NULL,
                        process = c("longitudinal", "event"),
                        type_pred = c("response", "link"),
                        type = c("subject_specific", "mean_subject"),
                        control = NULL, ...) {
    process <- match.arg(process)
    type_pred <- match.arg(type_pred)
    type <- match.arg(type)
    con <- list(all_times = FALSE, times_per_id = FALSE, level = 0.95,
                return_newdata = FALSE, use_Y = TRUE, return_mcmc = FALSE,
                n_samples = 200L, n_mcmc = 55L, parallel = "snow",
                cores = NULL, seed = 123L)
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0) {
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    }
    id_var <- object$model_info$var_names$idVar
    time_var <- object$model_info$var_names$time_var
    Time_var <- object$model_info$var_names$Time_var
    event_var <- object$model_info$var_names$event_var
    type_censoring <- object$model_info$type_censoring
    respVars <- unlist(object$model_info$var_names$respVars)
    check_varNames <- function (object, newdata, id_var,
                                process = c("Event", "Longitudinal")) {
        process <- match.arg(process)
        name_data <- deparse(substitute(newdata))
        if (process == "Event") {
            vars_S <- c(all.vars(object$model_info$terms$terms_Surv), id_var)
            missing_vars <- vars_S[!vars_S %in% names(newdata)]
            if (length(missing_vars)) {
                stop("the data.frame '", name_data, "' should contain the ",
                     "variable(s): ", paste(missing_vars, collapse = ", "),
                     ". \nThe '", paste(Time_var, collapse = ", "),
                     "' variable(s) should denote the last time the subjects",
                     " were event-free, and\nthe '", event_var,
                     "' variable should be set to 0.\n")
            }
        } else {
            termsL <- object$model_info$terms$terms_FE
            vars_L <- c(unlist(lapply(termsL, all.vars), use.names = FALSE), id_var)
            missing_vars <- vars_L[!vars_L %in% names(newdata)]
            if (length(missing_vars)) {
                stop("the data.frame '", name_data, "' should contain the ",
                     "variable(s): ", paste(missing_vars, collapse = ", "),
                     "\nFor the 'newdata2' or 'newdata2$newdataL', you will also ",
                     "need to provide the longitudinal\noutcome variables by ",
                     "setting them to a random value (these values are not used in ",
                     "the computations).\n")
            }
        }
    }
    if (object$model_info$CR_MS && is.data.frame(newdata)) {
        stop("for competing risks and multi-state models, argument 'newdata' ",
             "must be a list of two data.frames, one for the longitudinal ",
             "outcomes and one for the event process, the latter under the ",
             "correct format.\n")
    }
    if (!is.data.frame(newdata)) {
        if (!is.list(newdata) || length(newdata) != 2 ||
            !all(names(newdata) %in% c("newdataL", "newdataE"))) {
            stop("'newdata' must be a list with two data.frame elements ",
                 "named 'newdataL' and 'newdataE'.\n")
        }
        check_varNames(object, newdata$newdataE, id_var, "E")
        check_varNames(object, newdata$newdataL, id_var, "L")
        unq_ids_L <- newdata$newdataL[[id_var]]
        unq_ids_E <- newdata$newdataE[[id_var]]
        if (!all(unq_ids_L %in% unq_ids_E) || !all(unq_ids_E %in% unq_ids_L)) {
            stop("the subject id's in the datasets 'newdata$newdataL' and ",
                 "'newdata$newdataE' do not match.\n")
        }
    }
    if (is.data.frame(newdata)) {
        check_varNames(object, newdata, id_var, "E")
        check_varNames(object, newdata, id_var, "L")
    }
    if (!is.null(newdata2) && !is.data.frame(newdata2)) {
        if (!is.list(newdata2) || length(newdata2) != 2 ||
            !all(names(newdata2) %in% c("newdataL", "newdataE"))) {
            stop("'newdata2' must be a list with two data.frame elements ",
                 "named 'newdataL' and 'newdataE'.\n")
        }
        check_varNames(object, newdata2$newdataE, id_var, "E")
        check_varNames(object, newdata2$newdataL, id_var, "L")
        unq_ids_L <- newdata2$newdataL[[id_var]]
        unq_ids_E <- newdata2$newdataE[[id_var]]
        if (!all(unq_ids_L %in% unq_ids_E) || !all(unq_ids_E %in% unq_ids_L)) {
            stop("the subject id's in the datasets 'newdata2$newdataL' and ",
                 "'newdata2$newdataE' do not match.\n")
        }
    }
    if (!is.null(newdata2) && is.data.frame(newdata2)) {
        check_varNames(object, newdata2, id_var, "E")
        check_varNames(object, newdata2, id_var, "L")
    }
    if (is.null(con$cores)) {
        n <- if (!is.data.frame(newdata)) length(unique(newdata$newdataL[[id_var]]))
        else length(unique(newdata[[id_var]]))
        con$cores <- if (n > 20) 4L else 1L
    }
    components_newdata <-
        get_components_newdata(object, newdata, con$n_samples,
                               con$n_mcmc, con$parallel, con$cores, con$seed,
                               con$use_Y)
    if (process == "longitudinal") {
        predict_Long(object, components_newdata, newdata, newdata2, times,
                     con$all_times, con$times_per_id, type, type_pred, con$level,
                     con$return_newdata, con$return_mcmc)
    } else {
        predict_Event(object, components_newdata, newdata, newdata2, times,
                      con$times_per_id, con$level, con$return_newdata,
                      con$return_mcmc)
    }
}

plot.predict_jm <- function (x, x2 = NULL, subject = 1, outcomes = 1,
                             fun_long = NULL, fun_event = NULL,
                             CI_long = TRUE, CI_event = TRUE,
                             xlab = "Follow-up Time", ylab_long = NULL,
                             ylab_event = "Cumulative Risk", main = "",
                             lwd_long = 2, lwd_event = 2, ylim_event = c(0, 1),
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
    Time_var <- attr(x, "Time_var")
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
    if (!is.null(pred_Event)) xlim <- range(xlim, pred_Event[[Time_var]])
    plot_long_i <- function (outcome, add_xlab = FALSE, box = TRUE,
                             cex_axis = cex_axis) {
        ind <- pos_outcomes[outcome]
        outcome_i <- match(outcome, outcomes)
        f <- fun_long[[outcome_i]]
        preds <- f(pred_Long[[ind]])
        low <- f(pred_Long[[ind + 1]])
        upp <- f(pred_Long[[ind + 2]])
        times <- pred_Long[[time_var]]
        na_preds <- is.na(preds)
        preds <- preds[!na_preds]
        low <- low[!na_preds]
        upp <- upp[!na_preds]
        times <- times[!na_preds]
        ry <- range(preds, low, upp, na.rm = TRUE)
        ry <- range(ry[1L] * 0.8, ry[2L] * 1.2) # <---
        rx <- range(times, na.rm = TRUE)
        y_lim <- if (ylim_long_outcome_range) {
            range(f(ranges[[outcome]]), ry, na.rm = TRUE)
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
        unq_strata <- sort(unique(strata))
        col_line_event <- rep(col_line_event, length.out = length(unq_strata))
        fill_CI_event <- rep(fill_CI_event, length.out = length(unq_strata))
        times <- pred_Event[[Time_var]]
        ry <- sort(fun_event(c(0, 1)))
        rx <- range(times, na.rm = TRUE)
        plot(rx, ry, type = "n", xlab = "", ylab = "", xlim = xlim,
             axes = FALSE, col.axis = col_axis, col.lab = col_axis,
             ylim = ylim_event)
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
            pp <- par("usr")[3] + pos_ylab_long * diff(par("usr")[3:4])
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
            pp <- par("usr")[3] + pos_ylab_long * diff(par("usr")[3:4])
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

rc_setup <- function(rc_data, trm_data,
                     idVar = "id", statusVar = "status",
                     startVar = "start", stopVar = "stop",
                     trm_censLevel,
                     nameStrata = "strata", nameStatus = "status") {
  # warnings
  rc_bol <- c(idVar, statusVar, startVar, stopVar) %in% names(rc_data)
  if(any(!rc_bol)) {
    stop(paste0("\nThe variable '", c(idVar, statusVar, startVar, stopVar)[!rc_bol],
                "' is not present in 'rc_data' dataset."))
  }
  trm_bol <- c(idVar, statusVar, stopVar) %in% names(trm_data)
  if(any(!trm_bol)) {
    stop(paste0("\nThe variable '", c(idVar, statusVar, stopVar)[!trm_bol],
                "' is not present in 'trm_data' dataset."))
  }
  if(!setequal(rc_data[[idVar]],  trm_data[[idVar]])) {
    stop("The groups/subjects in both datasets do not seem to match.")
  }
  if(any(rc_data[[startVar]] > rc_data[[stopVar]])) {
    stop(paste0("'", stopVar, "' cannot be smaller than '", startVar," in the recurring event data.'"))
  }
  rc_data <- rc_data[order(rc_data[[idVar]], rc_data[[startVar]]), ]
  trm_data <- trm_data[order(trm_data[[idVar]]), ]
  if(any(rc_data[[stopVar]] > trm_data[[stopVar]][rc_data[[idVar]]])) {
    stop(paste0("'", stopVar, "' in the recurring event data cannot be larger than '", stopVar," in the terminal event data.'"))
  }
  # create new dataset
  ## CR dataset
  n <- nrow(trm_data)
  unqLevs <- unique(trm_data[[statusVar]])
  unqLevs <- unqLevs[unqLevs != trm_censLevel]
  status <- trm_data[[statusVar]] != trm_censLevel
  dataOut1 <- trm_data[rep(seq_len(n), each = length(unqLevs)), , drop = FALSE]
  dataOut1[[nameStrata]] <- rep(unqLevs, times = n)
  dataOut1[[nameStatus]] <- as.numeric(dataOut1[[statusVar]] == dataOut1[[nameStrata]])
  dataOut1[[startVar]] <- 0
  dataOut1[[nameStrata]] <- paste0("T", dataOut1[[nameStrata]])
  ## Rec dataset
  dataOut2 <- rc_data
  dataOut2[[nameStrata]] <- "R"
  dataOut2[[nameStatus]] <- dataOut2[[statusVar]]
  ## combine the 2 datasets
  dataOut <- rbind(dataOut1, dataOut2)
  dataOut[[nameStrata]] <- as.factor(dataOut[[nameStrata]]) # automatically assigns "R" as reference level based on the alphabetical order of the levels
  dataOut <- dataOut[order(dataOut[[idVar]], dataOut[[nameStrata]], dataOut[[startVar]]), ]
  rownames(dataOut) <- seq_len(nrow(dataOut))
  dataOut
}

predict.jmList <- function (object, weights, newdata = NULL, newdata2 = NULL,
                        times = NULL, process = c("longitudinal", "event"),
                        type_pred = c("response", "link"),
                        type = c("subject_specific", "mean_subject"),
                        control = NULL, ...) {
    process <- match.arg(process)
    type_pred <- match.arg(type_pred)
    type <- match.arg(type)
    con <- list(all_times = FALSE, times_per_id = FALSE, level = 0.95,
                return_newdata = FALSE, use_Y = TRUE, return_mcmc = FALSE,
                n_samples = 200L, n_mcmc = 55L, parallel = "snow",
                cores = parallelly::availableCores(omit = 1L), seed = 123L)
    control <- c(control, list(...))
    namC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% namC]) > 0) {
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    }
    obj <- object[[1L]]
    id_var <- obj$model_info$var_names$idVar
    time_var <- obj$model_info$var_names$time_var
    Time_var <- obj$model_info$var_names$Time_var
    event_var <- obj$model_info$var_names$event_var
    type_censoring <- obj$model_info$type_censoring
    respVars <- unlist(obj$model_info$var_names$respVars)
    if (obj$model_info$CR_MS && is.data.frame(newdata)) {
        stop("for competing risks and multi-state models, argument 'newdata' ",
             "must be a list of two data.frames, one for the longitudinal ",
             "outcomes and one for the event process, the latter under the ",
             "correct format.\n")
    }
    if (!is.data.frame(newdata)) {
        if (!is.list(newdata) || length(newdata) != 2 ||
            !all(names(newdata) %in% c("newdataL", "newdataE"))) {
            stop("'newdata' must be a list with two data.frame elements ",
                 "named 'newdataL' and 'newdataE'.\n")
        }
        for (i in seq_along(respVars)) {
            v <- respVars[i]
            if (is.null(newdata$newdataE[[v]])) {
                newdata$newdataE[[v]] <- rep(0.1, nrow(newdata$newdataE))
            }
        }
        termsL <- obj$model_info$terms$terms_FE_noResp
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
        if (length(Time_var) > 1L) {
            if (is.null(newdata[[Time_var[1L]]])) {
                newdata[[Time_var[1L]]] <- 0
            }
            if (is.null(newdata[[Time_var[2L]]])) {
                last_time <- function (x) max(x, na.rm = TRUE) + 1e-06
                f <- factor(newdata[[id_var]], unique(newdata[[id_var]]))
                newdata[[Time_var[2L]]] <- ave(newdata[[time_var]], f,
                                               FUN = last_time)
            }
        } else {
            if (is.null(newdata[[Time_var]])) {
                last_time <- function (x) max(x, na.rm = TRUE) + 1e-06
                f <- factor(newdata[[id_var]], unique(newdata[[id_var]]))
                newdata[[Time_var]] <- ave(newdata[[time_var]], f, FUN = last_time)
            }
        }
        termsL <- obj$model_info$terms$terms_FE_noResp
        all_vars <- unlist(lapply(termsL, all.vars), use.names = FALSE)
        all_vars <- all_vars[!all_vars %in% time_var]
        all_vars <- c(all_vars, all.vars(obj$model_info$terms$terms_Surv_noResp))
        missing_vars <- all_vars[!all_vars %in% names(newdata)]
        if (length(missing_vars)) {
            stop("the data.frame 'newdata' should contain the ",
                 "variable(s): ", paste(missing_vars, collapse = ", "), ".\n")
        }
    }
    if (!is.null(newdata2) && !is.data.frame(newdata2)) {
        if (!is.list(newdata2) || length(newdata2) != 2 ||
            !all(names(newdata2) %in% c("newdataL", "newdataE"))) {
            stop("'newdata2' must be a list with two data.frame elements ",
                 "named 'newdataL' and 'newdataE'.\n")
        }
        for (i in seq_along(respVars)) {
            v <- respVars[i]
            if (is.null(newdata2$newdataE[[v]])) {
                newdata2$newdataE[[v]] <- rep(0.1, nrow(newdata2$newdataE))
            }
        }
        termsL <- obj$model_info$terms$terms_FE_noResp
        all_vars <- unlist(lapply(termsL, all.vars), use.names = FALSE)
        all_vars <- all_vars[!all_vars %in% time_var]
        missing_vars <- all_vars[!all_vars %in% names(newdata2$newdataE)]
        if (length(missing_vars)) {
            stop("the data.frame 'newdata2$newdataE' should contain the ",
                 "variable(s): ", paste(missing_vars, collapse = ", "), ".\n")
        }
        missing_vars <- all_vars[!all_vars %in% names(newdata2$newdataL)]
        if (length(missing_vars)) {
            stop("the data.frame 'newdata2$newdataL' should contain the ",
                 "variable(s): ", paste(missing_vars, collapse = ", "), ".\n")
        }
    }
    if (!is.null(newdata2) && is.data.frame(newdata2)) {
        if (is.null(newdata2[[event_var]])) newdata2[[event_var]] <- 0
        if (length(Time_var) > 1L) {
            if (is.null(newdata2[[Time_var[1L]]])) {
                newdata2[[Time_var[1L]]] <- 0
            }
            if (is.null(newdata2[[Time_var[2L]]])) {
                last_time <- function (x) max(x, na.rm = TRUE) + 1e-06
                f <- factor(newdata2[[id_var]], unique(newdata2[[id_var]]))
                newdata2[[Time_var[2L]]] <- ave(newdata2[[time_var]], f,
                                                FUN = last_time)
            }
        } else {
            if (is.null(newdata2[[Time_var]])) {
                last_time <- function (x) max(x, na.rm = TRUE) + 1e-06
                f <- factor(newdata2[[id_var]], unique(newdata2[[id_var]]))
                newdata2[[Time_var]] <- ave(newdata2[[time_var]], f, FUN = last_time)
            }
        }
        termsL <- obj$model_info$terms$terms_FE_noResp
        all_vars <- unlist(lapply(termsL, all.vars), use.names = FALSE)
        all_vars <- all_vars[!all_vars %in% time_var]
        all_vars <- c(all_vars, all.vars(obj$model_info$terms$terms_Surv_noResp))
        missing_vars <- all_vars[!all_vars %in% names(newdata2)]
        if (length(missing_vars)) {
            stop("the data.frame 'newdata2' should contain the ",
                 "variable(s): ", paste(missing_vars, collapse = ", "), ".\n")
        }
    }
    cores <- min(con$cores, length(object))
    if (cores > 1L) {
        have_mc <- have_snow <- FALSE
        if (con$parallel == "multicore") {
            have_mc <- .Platform$OS.type != "windows"
        } else if (con$parallel == "snow") {
            have_snow <- TRUE
        }
        if (!have_mc && !have_snow) cores <- 1L
        loadNamespace("parallel")
    }
    if (cores > 1L) {
        if (have_mc) {
            preds <-
                parallel::mclapply(object, predict, newdata = newdata,
                                   newdata2 = newdata2, times = times,
                                   all_times = con$all_times,
                                   times_per_id = con$times_per_id,
                                   process = process, type_pred = type_pred,
                                   type = type, level = con$level,
                                   n_samples = con$n_samples,
                                   n_mcmc = con$n_mcmc,
                                   return_newdata = con$return_newdata,
                                   return_mcmc = TRUE, mc.cores = cores)
        } else {
            cl <- parallel::makePSOCKcluster(rep("localhost", cores))
            invisible(parallel::clusterEvalQ(cl, library("JMbayes2")))
            preds <-
                parallel::parLapply(cl, object, predict, newdata = newdata,
                                    newdata2 = newdata2, times = times,
                                    all_times = con$all_times,
                                    times_per_id = con$times_per_id,
                                    process = process, type_pred = type_pred,
                                    type = con$type, level = con$level,
                                    n_samples = con$n_samples,
                                    n_mcmc = con$n_mcmc,
                                    return_newdata = con$return_newdata,
                                    return_mcmc = TRUE)
            parallel::stopCluster(cl)
        }
    } else {
        preds <-
            lapply(object, predict, newdata = newdata,
                   newdata2 = newdata2, times = times,
                   all_times = con$all_times, times_per_id = con$times_per_id,
                   process = process, type_pred = type_pred,
                   type = type, level = con$level, n_samples = con$n_samples,
                   n_mcmc = con$n_mcmc, return_newdata = con$return_newdata,
                   return_mcmc = TRUE)
    }
    extract_mcmc <- function (x) {
        if (is.data.frame(x)) attr(x, "mcmc") else x[["mcmc"]]
    }
    MCMC <- lapply(preds, extract_mcmc)
    alp <- 1 - con$level
    if (is.list(MCMC[[1L]])) {
        n_outcomes <- length(MCMC[[1L]])
        pred_ <- qs <- vector("list", n_outcomes)
        for (j in seq_len(n_outcomes)) {
            weighted_MCMC <- Reduce("+", mapply2("*", MCMC[[j]], weights))
            pred_[[j]] <- rowMeans(weighted_MCMC)
            qs[[j]] <- matrixStats::rowQuantiles(weighted_MCMC,
                                                 probs = c(alp/2, 1 - alp/2))
        }
        names(pred_) <- names(qs) <- names(MCMC[[1L]])
        low <- lapply(qs, function (x) x[, 1L])
        upp <- lapply(qs, function (x) x[, 2L])
    } else {
        weighted_MCMC <- Reduce("+", mapply2("*", MCMC, weights))
        pred_ <- rowMeans(weighted_MCMC)
        qs <- matrixStats::rowQuantiles(weighted_MCMC,
                                        probs = c(alp/2, 1 - alp/2))
    }
    out <- preds[[1]]
    if (process == "event") {
        if (!is.data.frame(out)) {
            out$pred <- pred_
            out$low <- qs[, 1L]
            out$upp <- qs[, 2L]
        } else {
            out$pred_CIF <- pred_
            out$low_CIF <- qs[, 1L]
            out$upp_CIF <- qs[, 2L]
        }
    } else {
        if (!is.data.frame(out)) {
            out$preds <- pred_
            out$low <- low
            out$upp <- upp
        } else {
            ind <- grep("pred_", names(out), fixed = TRUE)
            out[ind] <- pred_
            ind <- grep("low_", names(out), fixed = TRUE)
            out[ind] <- low
            ind <- grep("upp_", names(out), fixed = TRUE)
            out[ind] <- upp
        }
    }
    out
}

simulate.jm <- function (object, nsim = 1L, seed = NULL, newdata = NULL,
                         process = c("longitudinal", "event"),
                         random_effects = c("posterior_means", "mcmc", "prior"),
                         Fforms_fun = NULL, include_outcome = FALSE,
                         tol = 0.001, iter = 100L, ...) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        runif(1L)
    if (is.null(seed))
        RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
        R.seed <- get(".Random.seed", envir = .GlobalEnv)
        set.seed(seed)
        RNGstate <- structure(seed, kind = as.list(RNGkind()))
        on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
    process <- match.arg(process)
    if (is.character(random_effects)) {
        random_effects <- match.arg(random_effects)
    }
    # information from fitted joint model
    ind_RE <- object$model_data$ind_RE
    control <- object$control
    if (is.null(newdata)) {
        n <- object$model_data$n
        y <- object$model_data$y
        X <- object$model_data$X
        Z <- object$model_data$Z
        idL_lp <- object$model_data$idL_lp
        times <- object$model_data$dataL[[object$model_info$var_names$time_var]]
    } else {
        id_var <- object$model_info$var_names$idVar
        if (is.null(id_var)) {
            stop("The id variable '", id_var, "' cannot be found in newdata.\n")
        }
        id <- newdata[[id_var]]
        id <- factor(id, levels = unique(id))
        n <- length(unique(id))
        terms_FE <- terms(object)
        frames_FE <- lapply(terms_FE, model.frame.default, data = newdata)
        terms_RE <- terms(object, type = "random")
        frames_RE <- lapply(terms_RE, model.frame.default, data = newdata)
        NAs_FE <- lapply(frames_FE, attr, "na.action")
        NAs_RE <- lapply(frames_RE, attr, "na.action")
        frames_FE <- mapply2(fix_NAs_fixed, frames_FE, NAs_FE, NAs_RE)
        frames_RE <- mapply2(fix_NAs_random, frames_RE, NAs_RE, NAs_FE)
        y <- lapply(frames_FE, model.response)
        y <- lapply(y, function (yy) {
            if (is.factor(yy)) as.numeric(yy != levels(yy)[1L]) else yy
        })
        X <- lapply(frames_FE, model.matrix.default, data = newdata)
        Z <- lapply(frames_RE, model.matrix.default, data = newdata)
        unq_id <- unique(id)
        id <- mapply2(exclude_NAs, NAs_FE, NAs_RE, MoreArgs = list(id = id))
        id[] <- lapply(id, match, table = unq_id)
        idL_lp <- lapply(id, function (x) match(x, unique(x)))
        times <- newdata[[object$model_info$var_names$time_var]]
    }
    n_outcomes <- length(idL_lp)
    families <- object$model_info$families
    has_sigmas <- as.logical(object$model_data$has_sigmas)
    Bspline_dgr <- control$Bsplines_degree
    knots <- control$knots
    basis <- control$basis
    timescale_base_hazard <- control$timescale_base_hazard
    if (is.null(newdata)) {
        dataS <- object$model_data$dataS
        Times <- object$model_data$Time_right
        event <- object$model_data$delta
        W <- object$model_data$W_h
        strt <- object$model_data$strata
    } else {
        dataS <- newdata[tapply(row.names(newdata), id, tail, n = 1L), ]
        terms_event <- terms(object, process = "event")
        frames_event <- model.frame.default(terms_event, data = dataS)
        Surv_Response <- model.response(frames_event)
        type_censoring <- attr(Surv_Response, "type")
        if (type_censoring == "right") {
            Times <- unname(Surv_Response[, "time"])
            event <- unname(Surv_Response[, "status"])
        } else if (type_censoring == "counting") {
            Times <- unname(Surv_Response[, "stop"])
            event <-  unname(Surv_Response[, "status"])
        } else {
            stop("simulate.jm() does not yet work for this type of censoring.\n")
        }
        W <- model.matrix.default(frames_event, data = dataS)[, -1L, drop = FALSE]
        ind_strata <- attr(terms_event, "specials")$strata
        strt <- if (is.null(ind_strata)) {
            rep(1, nrow(frames_event))
        } else {
            unclass(frames_event[[ind_strata]])
        }
    }
    # MCMC results
    ncz <- sum(sapply(Z, ncol))
    ind_betas <- grep("betas", names(object$statistics$Mean), fixed = TRUE)
    mcmc_betas <- object$mcmc[ind_betas]
    mcmc_betas[] <- lapply(mcmc_betas, function (x) do.call('rbind', x))
    mcmc_sigmas <- matrix(0.0, nrow(mcmc_betas[[1]]), n_outcomes)
    if (any(has_sigmas)) mcmc_sigmas[, has_sigmas] <- do.call('rbind', object$mcmc$sigmas)
    mcmc_bs_gammas <- do.call('rbind', object$mcmc$bs_gammas)
    has_gammas <- !is.null(object$mcmc$gammas)
    if (has_gammas) mcmc_gammas <- do.call('rbind', object$mcmc$gammas)
    mcmc_W_std_gammas <- do.call('rbind', object$mcmc$W_std_gammas)
    mcmc_alphas <- do.call('rbind', object$mcmc$alphas)
    mcmc_Wlong_std_alphas <- do.call('rbind', object$mcmc$Wlong_std_alphas)
    # random effects
    b <- ranef(object)
    if (length(random_effects) == 1L && random_effects == "mcmc") {
        mcmc_RE <- dim(object$mcmc[["b"]][[1L]])[3L] > 1L
        if (mcmc_RE) {
            mcmc_b <- abind::abind(object$mcmc[["b"]])
        } else {
            stop("refit the model using 'jm(..., save_random_effects = TRUE)'.\n")
        }
    }
    get_D <- function (x, n) {
        m <- matrix(0.0, n, n)
        m[lower.tri(m, TRUE)] <- x
        m <- m + t(m)
        diag(m) <- diag(m) * 0.5
        m
    }
    xx <- do.call('rbind', object$mcmc$D)
    mcmc_D <- array(0.0, c(ncz, ncz, nrow(xx)))
    for (m in seq_len(nrow(mcmc_betas[[1L]]))) {
        mcmc_D[, , m] <- get_D(xx[m, ], ncz)
    }
    # simulate outcome vectors
    if (process == "longitudinal") {
        sim_fun <- function (family, n, mu, phi) {
            switch(family$family,
                   "gaussian" = rnorm(n, mu, phi),
                   "Student's-t" = mu + phi * rt(n, df = .df),
                   "binomial" = rbinom(n, .N, mu),
                   "poisson" = rpois(n, mu),
                   "negative binomial" = rnbinom(n, size = phi, mu = mu),
                   "beta" = rbeta(n, shape1 = mu * phi, shape2 = phi * (1.0 - mu)))
        }
        indices <- sample(nrow(mcmc_betas[[1]]), nsim)
        val <- vector("list", n_outcomes)
        names(val) <- object$model_info$var_names$respVars_form
        if (length(random_effects) == 1L && random_effects == "prior") {
            simulated_RE <- array(0.0, c(n, ncz, nsim))
            for (j in seq_len(nsim)) {
                simulated_RE[, , j] <-
                    MASS::mvrnorm(n, rep(0, ncz), mcmc_D[, , indices[j]])
            }
        }
        for (i in seq_len(n_outcomes)) {
            rep_y <- matrix(0.0, nrow(X[[i]]), nsim)
            for (j in seq_len(nsim)) {
                # parameters
                jj <- indices[j]
                betas <- lapply(mcmc_betas, function (x) x[jj, ])
                sigmas <- mcmc_sigmas[jj, ]
                bb <- if (length(random_effects) == 1L && is.character(random_effects)) {
                    switch(random_effects,
                           "mcmc" = mcmc_b[, , jj],
                           "posterior_means" = b,
                           "prior" = simulated_RE[, , j])
                } else {
                    random_effects
                }
                FE <- c(X[[i]] %*% betas[[i]])
                RE <- rowSums(Z[[i]] * bb[idL_lp[[i]], ind_RE[[i]]])
                eta <- FE + RE
                mu <- families[[i]]$linkinv(eta)
                if (families[[i]][['family']] == "binomial") {
                    .N <- if (NCOL(y[[i]]) == 1) rep(1, length(mu)) else y[[i]][, 2]
                }
                if (families[[i]][['family']] == "Student's-t") {
                    .df <- families[[i]][['df']]
                }
                rep_y[, j] <- sim_fun(families[[i]], length(mu), mu, sigmas[i])
            }
            colnames(rep_y) <- paste0("sim_", seq_len(nsim))
            val[[i]] <- rep_y
        }
        if (length(random_effects) == 1L && random_effects == "prior") {
            outT <- simulate(object, nsim = nsim, newdata = newdata,
                         process = "event", seed = seed, tol = tol, iter = iter,
                         random_effects = simulated_RE, Fforms_fun = Fforms_fun)
            for (i in seq_len(n_outcomes)) {
                for (j in seq_len(ncol(val[[i]]))) {
                    val[[i]][times > outT$Times[idL_lp[[i]], j], j] <- NA_real_
                }
            }
        }
        if (include_outcome) val <- c(val, list(outcome = y))
    } else {
        if (is.null(Fforms_fun) || !is.function(Fforms_fun)) {
            stop("you need to provide the 'Fforms_fun' function; ",
                 "see the examples in `?simulate.jm` for more information.\n")
        }
        if (length(unique(object$model_data$strata)) > 1L) {
            stop("'simulate.jm()' does not currently support stratified joint models.\n")
        }
        log_hazard <- function (time, subj) {
            if (length(time) != length(subj)) subjj <- rep(subj, length(time))
            tt <- if (timescale_base_hazard == "identity") time else log(time)
            W0 <- create_W0(tt, knots, Bspline_dgr + 1L, strt[subj], basis)
            log_h0 <- c(W0 %*% bs_gammas) - rescale_factor
            covariates <- if (has_gammas) c(W[subj, , drop = FALSE] %*% gammas) else 0.0
            long <- c(Fforms_fun(time, betas, bb[subjj, , drop = FALSE],
                                 dataS[subjj, ]) %*% alphas)
            log_h0 + covariates + long
        }
        invS <- function (time, subj, .K) {
            GQ <- gaussLegendre(.K)
            wk <- GQ$wk
            sk <- GQ$sk
            K <- length(sk)
            P <- time / 2
            st <- outer(P, sk + 1)
            id_GQ <- rep(subj, each = K)
            log_Pwk <- rep(log(wk), length(P)) + rep(log(P), each = K)
            log_integrand <- log_Pwk + log_hazard(c(t(st)), id_GQ)
            c(rowsum(exp(log_integrand), id_GQ, reorder = FALSE)) + log_u[subj]
        }
        nr_root <- function (interval, fn, gr, n, tol = 0.001, iter = 100L) {
            subjs <- seq_len(n)
            Low <- interval[1L]
            Upp <- interval[2L]
            low <- rep(Low, n)
            upp <- rep(Upp, n)
            fn_Low <- fn(low, subjs, 32)
            fn_Upp <- fn(upp, subjs, 32)
            simulated_times <- numeric(n)
            negUpp <- fn_Upp < 0
            simulated_times[negUpp] <- Upp
            converged <- vector("logical", n)
            converged[negUpp] <- TRUE
            subjs <- subjs[!negUpp]
            low <- low[!negUpp]
            upp <- upp[!negUpp]
            fn_Low <- fn_Low[!negUpp]
            fn_Upp <- fn_Upp[!negUpp]
            tt <- tt_old <- rep((Low + Upp) / 2, length(subjs))
            ffn <- fn(tt, subjs, 18)
            ggr <- exp(gr(tt, subjs))
            for (i in seq_len(iter)) {
                # check convergence
                check <- abs(ffn) < tol
                simulated_times[subjs[check]] <- tt[check]
                converged[subjs[check]] <- TRUE
                if (all(converged)) break
                subjs <- subjs[!check]
                tt <- tt[!check]
                tt_old <- tt_old[!check]
                low <- low[!check]
                upp <- upp[!check]
                ffn <- ffn[!check]
                ggr <- ggr[!check]
                fn_Low <- fn_Low[!check]
                fn_Upp <- fn_Upp[!check]
                # propose new value using Newton-Raphson
                tt <- tt - ffn / ggr
                out_of_int <- tt < Low | tt > Upp
                if (any(out_of_int)) {
                    tt[out_of_int] <- (low[out_of_int] + upp[out_of_int]) / 2
                }
                ffn <- fn(tt, subjs, 18)
                ggr <- exp(gr(tt, subjs))
                ind1 <- ffn < 0 & ffn > fn_Low
                fn_Low[ind1] <- ffn[ind1]
                low[ind1] <- tt[ind1]
                ind2 <- ffn > 0 & ffn < fn_Upp
                fn_Upp[ind2] <- ffn[ind2]
                upp[ind2] <- tt[ind2]
                simulated_times[subjs] <- tt
            }
            simulated_times
        }
        valT <- eventT <- matrix(0.0, n, nsim)
        indices <- sample(nrow(mcmc_betas[[1]]), nsim)
        for (j in seq_len(nsim)) {
            jj <- indices[j]
            betas <- lapply(mcmc_betas, function (x) x[jj, ])
            bb <- if (length(random_effects) == 1L && is.character(random_effects)) {
                switch(random_effects,
                       "mcmc" = mcmc_b[, , jj],
                       "posterior_means" = b,
                       "prior" = MASS::mvrnorm(n, rep(0, ncz), mcmc_D[, , jj]))
            } else {
                random_effects[, , j]
            }

            bs_gammas <- mcmc_bs_gammas[jj, ]
            if (has_gammas) gammas <- mcmc_gammas[jj, ]
            alphas <- mcmc_alphas[jj, ]
            Wlong_std_alphas <- mcmc_Wlong_std_alphas[jj, ]
            W_std_gammas <- mcmc_W_std_gammas[jj, ]
            rescale_factor <- Wlong_std_alphas + W_std_gammas
            Up <- max(Times) * 1.05
            log_u <- log(runif(n))
            rep_Times <- nr_root(c(0, Up), invS, log_hazard, n,
                                 tol = tol, iter = iter)
            valT[, j] <- rep_Times
            eventT[, j] <- as.numeric(rep_Times < Up)
        }
        colnames(valT) <- colnames(eventT) <- paste0("sim_", seq_len(nsim))
        val <- list(Times = valT, event = eventT)
        if (include_outcome) {
            val <- c(val, outcome = list("Times" = Times, "event" = event))
        }
    }
    attr(val, "seed") <- RNGstate
    val
}

ppcheck <- function (object, nsim = 40L, newdata = NULL, seed = 123L,
                     process = c("longitudinal", "event"),
                     outcomes = Inf, percentiles = c(0.025, 0.975),
                     random_effects = c("posterior_means", "mcmc", "prior"),
                     Fforms_fun = NULL, ...) {
    process <- match.arg(process)
    random_effects <- match.arg(random_effects)
    trapezoid_rule <- function (f, x) {
        sum(0.5 * diff(x) * (f[-length(x)] + f[-1L]))
    }
    bind <- function (sims) {
        out <- sims[[1L]]
        for (i in seq_along(sims)[-1L]) {
            out_i <- sims[[i]]
            for (j in seq_along(out)) {
                out[[j]] <- if (is.matrix(out[[j]]))
                    rbind(out[[j]], out_i[[j]]) else c(out[[j]], out_i[[j]])
            }
        }
        out
    }
    if (process == "longitudinal") {
        out <- if (inherits(object, 'list') && inherits(object[[1L]], 'jm') &&
                   inherits(newdata, 'list')) {
            sims_per_fold <-
                mapply(simulate, object = object, newdata = newdata,
                       MoreArgs = list(process = "longitudinal",
                                       include_outcome = TRUE,
                                       random_effects = random_effects,
                                       seed = seed, nsim = nsim,
                                       Fforms_fun = Fforms_fun, ...),
                       SIMPLIFY = FALSE)
            bind(sims_per_fold)
        } else {
            simulate(object, nsim = nsim, newdata = newdata,
                     process = "longitudinal", include_outcome = TRUE,
                     random_effects = random_effects, seed = seed,
                     Fforms_fun = Fforms_fun, ...)
        }
        yy <- out$outcome
        out <- out[names(out) != "outcome"]
        n_outcomes <- length(out)
        index <- seq_len(n_outcomes)
        if (outcomes < Inf) index <- index[index %in% outcomes]
        for (j in index) {
            y <- yy[[j]]
            r1 <- quantile(y, probs = percentiles[1L], na.rm = TRUE)
            r2 <- quantile(y, probs = percentiles[2L], na.rm = TRUE)
            x_vals <- seq(r1, r2, length.out = 600)
            rep_y <- apply(out[[j]], 2L, function (x, x_vals) ecdf(x)(x_vals),
                            x_vals = x_vals)
            F0 <- ecdf(y)
            # Dvoretzky–Kiefer–Wolfowitz inequality
            # https://stats.stackexchange.com/questions/181724/confidence-intervals-for-ecdf
            #xx <- get("x", envir = environment(F0))
            #ff <- get("y", envir = environment(F0))
            #eps <- sqrt(log(2 / 0.05) / (2 * length(y)))
            #F0u <- pmin(ff + eps, 1)
            #F0l <- pmax(ff - eps, 0)
            #F0u <- stepfun(xx, c(F0u, 1))(x_vals)
            #F0l <- stepfun(xx, c(F0l, F0l[length(y)]))(x_vals)
            F0 <- F0(x_vals)
            se <- sqrt(F0 * (1 - F0) / length(y))
            F0u <- pmin(F0 + 1.959964 * se, 1)
            F0l <- pmax(F0 - 1.959964 * se, 0)
            MISE <- mean(apply((rep_y - F0)^2, 2L, trapezoid_rule, x = x_vals))
            matplot(x_vals, rep_y, type = "s", lty = 1, col = "lightgrey",
                    xlab = object$model_info$var_names$respVars_form[[j]],
                    ylab = "Empirical CDF", ylim = c(0, 1))
            lines(x_vals, F0, lwd = 1.5, type = "s")
            lines(x_vals, F0l, lwd = 1.5, lty = 2, type = "s")
            lines(x_vals, F0u, lwd = 1.5, lty = 2, type = "s")
            legend("bottomright", c("replicated data", "observed data"),
                   lty = 1, col = c("lightgrey", "black"), bty = "n", cex = 0.9)
            rootMISE <- round(sqrt(MISE), 5)
            text(r1 + 0.15 * (r2 - r1), 0.9, bquote(sqrt(MISE) == .(rootMISE)))
        }
    } else {
        out <- if (inherits(object, 'list') && inherits(object[[1L]], 'jm') &&
                   inherits(newdata, 'list')) {
            sims_per_fold <-
                mapply(simulate, object = object, newdata = newdata,
                       MoreArgs = list(process = "event", include_outcome = TRUE,
                                       Fforms_fun = Fforms_fun,
                                       random_effects = random_effects,
                                       seed = seed, nsim = nsim, ...),
                       SIMPLIFY = FALSE)
            bind(sims_per_fold)
        } else {
            simulate(object, nsim = nsim, newdata = newdata,
                     process = "event", seed = seed, include_outcome = TRUE,
                     random_effects = random_effects, Fforms_fun = Fforms_fun,
                     ...)
        }
        Times <- out$outcome.Times
        event <- out$outcome.event
        out <- out[names(out) != "outcome"]
        r2 <- quantile(Times, probs = percentiles[2L], na.rm = TRUE)
        x_vals <- seq(0, r2, length.out = 500)
        rep_T <- apply(out$Times, 2L, function (x, x_vals) ecdf(x)(x_vals),
                        x_vals = x_vals)
        ss <- summary(survfit(Surv(Times, event) ~ 1), times = x_vals)
        F0 <- 1 - ss$surv
        F0_low <- 1 - ss$lower
        F0_upp <- 1 - ss$upper
        MISE <- mean(apply((rep_T - F0)^2, 2L, trapezoid_rule, x = x_vals))
        matplot(x_vals, rep_T, type = "s", lty = 1, col = "lightgrey",
                xlab = "Times", ylab = "Empirical CDF", ylim = c(0, 1))
        lines(x_vals, F0, lwd = 1.5, type = "s")
        lines(x_vals, F0_low, lty = 2, lwd = 1.5, type = "s")
        lines(x_vals, F0_upp, lty = 2, lwd = 1.5, type = "s")
        legend("bottomright", c("replicated data", "observed data"),
               lty = 1, col = c("lightgrey", "black"), bty = "n", cex = 0.9)
        rootMISE <- round(sqrt(MISE), 5)
        text(0.15 * r2, 0.9, bquote(sqrt(MISE) == .(rootMISE)))
    }
}


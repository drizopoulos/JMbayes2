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

coef.summary.jm <- function (object, process = c("longitudinal", "event"), ...) {
    process <- match.arg(process)
    if (process == "event") {
        object[['Survival']]
    } else {
        ind <- grep("Outcome", names(object), fixed = TRUE)
        out <- object[ind]
        names(out) <- object[['respVars']]
        out
    }
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

fitted.jm <- function (object, ...) {
    fits <- function (X, betas, Z, b, id, ind_RE) {
        FE <- c(X %*% betas)
        RE <- rowSums(Z * b[id, ind_RE, drop = FALSE])
        FE + RE
    }
    means <- object$statistics$Mean
    ind <- grep("betas", names(means), fixed = TRUE)
    out <-
        mapply2(fits, X = object$model_data[['X']], betas = means[ind],
                Z = object$model_data[['Z']], id = object$model_data$idL_lp,
                ind_RE = object$model_data$ind_RE, MoreArgs = list(b = means[['b']]))
    names(out) <- unlist(object$model_info$var_names$respVars)
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
                           linewidth = 1, alpha = 0.8,
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
                    linewidth = linewidth, alpha = alpha) +
          ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
          theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
          scale_color_manual(values = ggcolthemes[[coltheme]]) +
          guides(color = guide_legend(override.aes = list(alpha = 1)))
      } else {
        gplots[[i]] <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
          geom_line(aes(x = iteration, y = value, color = chain),
                    linewidth = linewidth, alpha = alpha) +
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
                    linewidth = linewidth, alpha = alpha) +
          ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
          theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
          scale_color_manual(values = ggcolthemes[[coltheme]]) +
          guides(color = guide_legend(override.aes = list(alpha = 1)))
        print(g)
      } else {
        g <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
          geom_line(aes(x = iteration, y = value, color = chain),
                    linewidth = linewidth, alpha = alpha) +
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
                      linewidth = 1, alpha = 0.6,
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
                           linewidth = linewidth, alpha = alpha) +
              ggtitle(paste('Density plot of ', unique(ggdata$parm)[i])) +
              theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
              scale_color_manual(values = ggcolthemes[[coltheme]]) +
              scale_fill_manual(values = ggcolthemes[[coltheme]]) +
              guides(color = guide_legend(override.aes = list(alpha = 1)))
          } else {
            gplots[[i]] <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
              geom_density(aes(x = value, color = chain, fill = chain),
                           linewidth = linewidth, alpha = alpha) +
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
                           linewidth = linewidth, alpha = alpha) +
              ggtitle(paste('Density plot of ', unique(ggdata$parm)[i])) +
              theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
              scale_color_manual(values = ggcolthemes[[coltheme]]) +
              scale_fill_manual(values = ggcolthemes[[coltheme]]) +
              guides(color = guide_legend(override.aes = list(alpha = 1)))
            print(g)
          } else {
            g <- ggplot(ggdata[ggdata$parm %in% unique(ggdata$parm)[i], ]) +
              geom_density(aes(x = value, color = chain, fill = chain),
                           linewidth = linewidth, alpha = alpha) +
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
                return_params_mcmc = FALSE, n_samples = 200L, n_mcmc = 55L,
                parallel = "snow", cores = NULL, seed = 123L)
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
                     con$return_newdata, con$return_mcmc, con$return_params_mcmc)
    } else {
        predict_Event(object, components_newdata, newdata, newdata2, times,
                      con$times_per_id, con$level, con$return_newdata,
                      con$return_mcmc, con$return_params_mcmc)
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
  if (any(!rc_bol)) {
    stop(paste0("The variable '", c(idVar, statusVar, startVar, stopVar)[!rc_bol],
                "' is not present in 'rc_data' dataset."))
  }
  trm_bol <- c(idVar, statusVar, stopVar) %in% names(trm_data)
  if (any(!trm_bol)) {
    stop(paste0("The variable '", c(idVar, statusVar, stopVar)[!trm_bol],
                "' is not present in 'trm_data' dataset."))
  }
  if (!setequal(rc_data[[idVar]], trm_data[[idVar]])) {
    stop("The groups/subjects in both datasets do not seem to match.")
  }
  if (any(rc_data[[startVar]] > rc_data[[stopVar]])) {
    stop(paste0("'", stopVar, "' cannot be smaller than '", startVar,
                "' in the recurring event data."))
  }
  rc_data  <- rc_data[order(rc_data[[idVar]], rc_data[[startVar]]), ]
  trm_data <- trm_data[order(trm_data[[idVar]]), ]
  id_match <- match(rc_data[[idVar]], trm_data[[idVar]])
  if (any(rc_data[[stopVar]] > trm_data[[stopVar]][id_match])) {
    stop(paste0("'", stopVar, "' in the recurring event data cannot be larger than '",
                stopVar,"' in the terminal event data."))
  }
  # create new dataset
  ## CR dataset
  n <- nrow(trm_data)
  unqLevs <- unique(trm_data[[statusVar]])
  unqLevs <- unqLevs[unqLevs != trm_censLevel]
  dataOut1 <- trm_data[rep(seq_len(n), each = length(unqLevs)), , drop = FALSE]
  dataOut1[[nameStrata]] <- rep(unqLevs, times = n)
  dataOut1[[nameStatus]] <- as.numeric(dataOut1[[statusVar]] == dataOut1[[nameStrata]])
  dataOut1[[startVar]] <- 0
  dataOut1[[nameStrata]] <- paste0("T", dataOut1[[nameStrata]])
  ## Rec dataset
  dataOut2 <- rc_data
  dataOut2[[nameStrata]] <- "R"
  dataOut2[[nameStatus]] <- as.numeric(dataOut2[[statusVar]])
  ## combine the 2 datasets
  names1 <- names(dataOut1)
  names2 <- names(dataOut2)
  common_names <- intersect(names1, names2)
  miss1 <- setdiff(names2, names1)
  miss2 <- setdiff(names1, names2)
  if (length(miss1)) {
    warning("The following variables were missing in the 'trm_data' and were created as NA: ",
            paste(miss1, collapse = ", "), ".")
    dataOut1[miss1] <- NA # add missing vars as NA
  }
  if (length(miss2)) {
    warning("The following variables were missing in the 'rc_data' and were created as NA: ",
            paste(miss2, collapse = ", "), ".")
    dataOut2[miss2] <- NA
  }
  class1 <- sapply(dataOut1[common_names], class)
  class2 <- sapply(dataOut2[common_names], class)
  conflict_names <- common_names[class1 != class2]
  if (length(conflict_names)) {
    warning("The following variables had the same name but different classes and were renamed: ",
            paste(conflict_names, collapse = ", "), ".")
    names(dataOut1)[match(conflict_names, names(dataOut1))] <- paste0(conflict_names,
                                                                      "_trm") # rename vars that are in BOTH and have different class
    names(dataOut2)[match(conflict_names, names(dataOut2))] <- paste0(conflict_names,
                                                                      "_rc")
    dataOut1[paste0(conflict_names, "_rc")] <- NA
    dataOut2[paste0(conflict_names, "_trm")] <- NA
  }
  dataOut2 <- dataOut2[names(dataOut1)] # reorder columns in the same order
  dataOut <- rbind(dataOut1, dataOut2)
  dataOut[[nameStrata]] <- as.factor(dataOut[[nameStrata]]) # automatically assigns "R" as reference level based on the alphabetical order of the levels
  dataOut <- dataOut[order(dataOut[[idVar]],
                           dataOut[[nameStrata]],
                           dataOut[[startVar]]), ]
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
                         params_mcmc = NULL, Fforms_fun = NULL,
                         include_outcome = FALSE, return_random_effects = FALSE,
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
        times <- mapply2(exclude_NAs, object$model_info$NAs_FE,
                         object$model_info$NAs_RE, MoreArgs = list(id = times))
    } else {
        id_var <- object$model_info$var_names$idVar
        if (is.null(id_var)) {
            stop("The id variable '", id_var, "' cannot be found in newdata.\n")
        }
        id <- newdata[[id_var]]
        id <- id. <- factor(id, levels = unique(id))
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
        X <- mapply2(model.matrix.default, terms_FE, frames_FE)
        Z <- mapply2(model.matrix.default, terms_RE, frames_RE)
        unq_id <- unique(id)
        id <- mapply2(exclude_NAs, NAs_FE, NAs_RE, MoreArgs = list(id = id))
        id[] <- lapply(id, match, table = unq_id)
        idL_lp <- lapply(id, function (x) match(x, unique(x)))
        times <- newdata[[object$model_info$var_names$time_var]]
        times <- mapply2(exclude_NAs, NAs_FE, NAs_RE,
                         MoreArgs = list(id = times))
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
        dataS <- newdata[tapply(row.names(newdata), id., tail, n = 1L), ]
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
    if (is.null(params_mcmc)) {
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
    } else {
        mcmc_betas <- params_mcmc[['betas']]
        mcmc_sigmas <- params_mcmc[['sigmas']]
        mcmc_bs_gammas <- params_mcmc[['bs_gammas']]
        mcmc_gammas <- params_mcmc[['gammas']]
        mcmc_W_std_gammas <- params_mcmc[['W_std_gammas']]
        mcmc_alphas <- params_mcmc[['alphas']]
        mcmc_Wlong_std_alphas <- params_mcmc[['Wlong_std_alphas']]
        mcmc_b <- params_mcmc[['b']]
        mcmc_D <- params_mcmc[['D']]
        if (nsim > dim(mcmc_b)[3L]) {
            warning("'nsim' is greater than the available MCMC sample size. ",
                    "It is set to ", dim(mcmc_b)[3L], ".\n")
            nsim <- dim(mcmc_b)[3L]
        }
        if (length(random_effects) == 1L && is.character(random_effects) &&
            random_effects == "posterior_means") {
            stop("'random_effects' cannot be set to 'posterior_means' ",
                 "when 'params_mcmc' is not NULL.\nYou will need to include ",
                 "the MCMC sample of random effects in the 'params_mcmc' list.\n")
        }
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
                   "beta" = rbeta(n, shape1 = mu * phi, shape2 = phi * (1 - mu)),
                   "Gamma" = rgamma(n, shape = phi, scale = mu / phi))
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
            if (include_outcome) {
                attr(y[[i]], "times") <- times[[i]]
                attr(y[[i]], "id") <- idL_lp[[i]]
            }
            if (return_random_effects) {
                out_RE <- array(0.0, c(n, ncz, nsim))
            }
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
                    if (length(dim(random_effects)) > 2)
                        random_effects[, , j] else random_effects
                }
                if (return_random_effects) {
                    out_RE[, , j] <- bb
                }
                FE <- c(X[[i]] %*% betas[[i]])
                if (!is.matrix(bb)) bb <- rbind(bb)
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
            if (return_random_effects) {
                attr(val[[i]], "RE") <- out_RE
            }
        }
        if (length(random_effects) == 1L && random_effects == "prior") {
            outT <- simulate(object, nsim = nsim, newdata = newdata,
                         process = "event", seed = seed, tol = tol, iter = iter,
                         random_effects = simulated_RE, Fforms_fun = Fforms_fun)
            for (i in seq_len(n_outcomes)) {
                for (j in seq_len(ncol(val[[i]]))) {
                    na_ind <- times[[i]] > outT$Times[idL_lp[[i]], j]
                    val[[i]][na_ind, j] <- NA_real_
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
            subjj <- if (length(time) != length(subj)) rep(subj, length(time)) else subj
            tt <- if (timescale_base_hazard == "identity") time else log(time)
            W0 <- create_W0(tt, knots, Bspline_dgr, strt[subj], basis,
                            timescale_base_hazard)
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
                     process = c("longitudinal", "event", "joint"),
                     type = c("ecdf", "average-evolution", "variance-function",
                              "variogram", "surv-uniform"),
                     CI_ecdf = c("none", "binomial", "Dvoretzky-Kiefer-Wolfowitz"),
                     CI_loess = FALSE,
                     outcomes = Inf, percentiles = c(0.025, 0.975),
                     random_effects = c("posterior_means", "mcmc", "prior"),
                     params_mcmc = NULL, Fforms_fun = NULL, plot = TRUE,
                     add_legend = TRUE, pos_legend = c("bottomright", "right"),
                     main = "", xlab = NULL, ylab = NULL,
                     col_obs = "black", col_rep = "lightgrey", lty_obs = 1,
                     lty_rep = 1, lwd_obs = 1.5, lwd_rep = 1, line_main = NA,
                     cex.main = 1.2, ylim = NULL, ...) {
    process <- match.arg(process)
    type <- match.arg(type)
    CI_ecdf <- match.arg(CI_ecdf)
    trapezoid_rule <- function (f, x) {
        sum(0.5 * diff(x) * (f[-length(x)] + f[-1L]))
    }
    bind <- function (sims) {
        Cbind <- function (x1, x2) {
            if (is.matrix(x1)) rbind(x1, x2) else c(x1, x2)
        }
        out <- sims[[1L]]
        for (i in seq_along(sims)[-1L]) {
            out_i <- sims[[i]]
            for (j in seq_along(out)) {
                if (is.list(out[[j]])) {
                    for (k in seq_along(out[[j]])) {
                        tt <- c(attr(out[[j]][[k]], "times"),
                                attr(out_i[[j]][[k]], "times"))
                        ii <- c(attr(out[[j]][[k]], "id"),
                                attr(out_i[[j]][[k]], "id"))
                        out[[j]][[k]] <- Cbind(out[[j]][[k]], out_i[[j]][[k]])
                        attr(out[[j]][[k]], "times") <- tt
                        attr(out[[j]][[k]], "id") <- ii
                    }
                } else out[[j]] <- Cbind(out[[j]], out_i[[j]])
            }
        }
        out
    }
    if (process == "longitudinal") {
        list_of_jms <- inherits(object, 'list') && inherits(object[[1L]], 'jm') &&
            inherits(newdata, 'list')
        out <- if (list_of_jms) {
            sims_per_fold <- if (is.null(params_mcmc)) {
                mapply2(simulate, object = object, newdata = newdata,
                        MoreArgs = list(process = "longitudinal",
                                        include_outcome = TRUE,
                                        random_effects = random_effects,
                                        seed = seed, nsim = nsim,
                                        Fforms_fun = Fforms_fun))
            } else {
                mapply2(simulate, object = object, newdata = newdata,
                        params_mcmc = params_mcmc,
                        MoreArgs = list(process = "longitudinal",
                                        include_outcome = TRUE,
                                        random_effects = random_effects,
                                        seed = seed, nsim = nsim,
                                        Fforms_fun = Fforms_fun))
            }
            bind(sims_per_fold)
        } else {
            simulate(object, nsim = nsim, newdata = newdata,
                     process = "longitudinal", include_outcome = TRUE,
                     random_effects = random_effects, params_mcmc = params_mcmc,
                     seed = seed, Fforms_fun = Fforms_fun)
        }
        yy <- out$outcome
        out <- out[names(out) != "outcome"]
        n_outcomes <- length(out)
        index <- seq_len(n_outcomes)
        if (outcomes < Inf) index <- index[index %in% outcomes]
        nindex <- length(index)
        plot_values <- if (plot) vector("list", nindex)
        if (is.null(xlab)) {
            xlab <- if (type == "ecdf") {
                    if (list_of_jms) object[[1L]]$model_info$var_names$respVars_form
                else object$model_info$var_names$respVars_form
                } else
                switch(type, 'average-evolution' = ,
                       'variance-function' = 'Follow-up Time',
                       'variogram' = 'Time Lags')
        }
        if (is.null(ylab)) {
            ylab <- if (type == "ecdf") "Empirical CDF" else
                switch(type, 'average-evolution' = 'Average Evolution',
                       'variance-function' =  expression(sqrt(abs("Std. Residuals"))),
                       'variogram' = 'Half Squared Differences')
        }
        xlab <- rep(xlab, length.out = nindex)
        ylab <- rep(ylab, length.out = nindex)
        main <- rep(main, length.out = nindex)
        for (j in index) {
            jj <- match(j, unique(index))
            if (type == "ecdf") {
                y <- yy[[j]]
                r1 <- quantile(y, probs = percentiles[1L], na.rm = TRUE)
                r2 <- quantile(y, probs = percentiles[2L], na.rm = TRUE)
                x_vals <- seq(r1, r2, length.out = 600)
                rep_y <- apply(out[[j]], 2L, function (x, x_vals) ecdf(x)(x_vals),
                               x_vals = x_vals)
                F0 <- ecdf(y)
                F0 <- F0(x_vals)
                if (CI_ecdf == "binomial") {
                    se <- sqrt(F0 * (1 - F0) / length(y))
                    F0u <- pmin(F0 + 1.959964 * se, 1)
                    F0l <- pmax(F0 - 1.959964 * se, 0)
                } else {
                    # Dvoretzky–Kiefer–Wolfowitz inequality
                    # https://stats.stackexchange.com/questions/181724/confidence-intervals-for-ecdf
                    eps <- sqrt(log(2 / 0.05) / (2 * length(y)))
                    F0u <- pmin(F0 + eps, 1)
                    F0l <- pmax(F0 - eps, 0)
                }
                MISE <- mean(apply((rep_y - F0)^2, 2L, trapezoid_rule, x = x_vals))
                rootMISE <- round(sqrt(MISE), 5)
                if (plot) {
                    if (is.null(ylim)) {
                        yylim <- c(0, min(1, max(rep_y, F0u) + 0.1))
                    } else yylim <- ylim
                    matplot(x_vals, rep_y, type = "s", lty = lty_rep, lwd = lwd_rep,
                            col = col_rep, xlab = xlab[jj], ylab = ylab[jj],
                            ylim = yylim, ...)
                    title(main = main[jj], line = line_main, cex.main = cex.main)
                    lines(x_vals, F0, type = "s", lwd = lwd_obs, lty = lty_obs, col = col_obs)
                    if (CI_ecdf != "none") {
                        lines(x_vals, F0l, type = "s", lwd = lwd_obs, lty = 2, col = col_obs)
                        lines(x_vals, F0u, type = "s", lwd = lwd_obs, lty = 2, col = col_obs)
                    }
                    if (add_legend) {
                        if (!is.na(pos_legend[1L])) {
                            if (CI_ecdf != "none") {
                                legend(pos_legend[1L],
                                       legend = c("replicated data", "observed data", "95% CI"),
                                       lty = c(1, 1, 2), col = c(col_rep, col_obs, col_obs),
                                       bty = "n", cex = 0.9)
                            } else {
                                legend(pos_legend[1L],
                                       legend = c("replicated data", "observed data"),
                                       lty = c(1, 1), col = c(col_rep, col_obs),
                                       bty = "n", cex = 0.9)

                            }
                        }
                        if (!is.na(pos_legend[2L])) {
                            legend(pos_legend[2L], bty = "n",
                                   legend = bquote(sqrt(MISE) == .(rootMISE)))
                        }
                    }
                } else {
                    plot_values[[jj]] <-
                        list(x_vals = x_vals, obs = F0, obs_low = F0l,
                             obs_upp = F0u, rep = rep_y,
                             rootMISE = if (add_legend) rootMISE)
                }
            } else {
                loess.smooth2 <- function (x, y) {
                    loess.smooth(x, y, degree = 2, span = 0.75,
                                 family = "gaussian", evaluation = 200)
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
                y <- yy[[j]]
                X <- attr(y, "X")
                tt <- attr(y, "times")
                id <- attr(y, "id")
                DF <- gof_fun(y, tt, id, type)
                obs_loess <- loess.smooth2(DF[, 1L], DF[, 2L])
                if (CI_loess) {
                    loess_fit <-
                        loess(y ~ x, data.frame(x = DF[, 1L], y = DF[, 2L]),
                              defree = 2, span = 0.75, family = "gaussian")
                    preds <- predict(loess_fit, data.frame(x = obs_loess$x),
                                     se = TRUE)
                    low <- preds$fit + qt(0.025, preds$df) * preds$se.fit
                    upp <- preds$fit + qt(0.975, preds$df) * preds$se.fit
                }
                rep_loess <- matrix(0, length(obs_loess$y), ncol(out[[j]]))
                for (i in seq_len(ncol(out[[j]]))) {
                    not_na <- !is.na(out[[j]][, i])
                    DF <- gof_fun(out[[j]][not_na, i], tt[not_na], id[not_na], type)
                    loess_rep_i <- loess.smooth2(DF[, 1L], DF[, 2L])
                    rep_loess[, i] <- loess_rep_i$y
                }
                MISE <- mean(apply((rep_loess - obs_loess$y)^2, 2L,
                                   trapezoid_rule, x = obs_loess$x))
                rootMISE <- round(sqrt(MISE), 5)
                if (plot) {
                    if (is.null(ylim)) {
                        yylim <- range(obs_loess$y, rep_loess)
                    } else yylim <- ylim
                    if (CI_loess) yylim <- range(yylim, low, upp)
                    matplot(obs_loess$x, rep_loess, type = "l", col = col_rep,
                            lty = lty_rep, lwd = lwd_rep, ylim = yylim,
                            xlab = xlab[jj], ylab = ylab[jj], ...)
                    title(main = main[jj], line = line_main, cex.main = cex.main)
                    lines(obs_loess, lwd = lwd_obs, lty = lty_obs, col = col_obs)
                    if (CI_loess) {
                        lines(obs_loess$x, low, lwd = lwd_obs, lty = 2, col = col_obs)
                        lines(obs_loess$x, upp, lwd = lwd_obs, lty = 2, col = col_obs)
                    }
                    if (add_legend) {
                        if (!is.na(pos_legend[1L])) {
                            if (CI_loess) {
                                legend(pos_legend[1L],
                                       legend = c("replicated data", "observed data",
                                                  "95% CI"), lty = c(1, 1, 2),
                                       col = c(col_rep, col_obs, col_obs), bty = "n",
                                       cex = 0.9)
                            } else {
                                legend(pos_legend[1L],
                                       legend = c("replicated data", "observed data"),
                                       lty = c(1, 1), col = c(col_rep, col_obs),
                                       bty = "n", cex = 0.9)
                            }
                        }
                        if (!is.na(pos_legend[2L])) {
                            legend(pos_legend[2L], bty = "n",
                                   legend = bquote(sqrt(MISE) == .(rootMISE)))
                        }
                    }
                } else {
                    plot_values[[jj]] <-
                        list(x_vals = obs_loess$x, obs = obs_loess$y,
                             obs_low = if (CI_loess) low,
                             obs_upp = if (CI_loess) upp, rep = rep_loess,
                             rootMISE = if (add_legend) rootMISE)
                }
            }
        }
    } else if (process == "event") {
        out <- if (inherits(object, 'list') && inherits(object[[1L]], 'jm') &&
                   inherits(newdata, 'list')) {
            sims_per_fold <- if (is.null(params_mcmc)) {
                mapply2(simulate, object = object, newdata = newdata,
                        MoreArgs = list(process = "event",
                                        include_outcome = TRUE,
                                        Fforms_fun = Fforms_fun,
                                        random_effects = random_effects,
                                        seed = seed, nsim = nsim))
            } else {
                mapply2(simulate, object = object, newdata = newdata,
                        params_mcmc = params_mcmc,
                        MoreArgs = list(process = "event",
                                        include_outcome = TRUE,
                                        Fforms_fun = Fforms_fun,
                                        random_effects = random_effects,
                                        seed = seed, nsim = nsim))
            }
            bind(sims_per_fold)
        } else {
            simulate(object, nsim = nsim, newdata = newdata,
                     process = "event", seed = seed, include_outcome = TRUE,
                     random_effects = random_effects, params_mcmc = params_mcmc,
                     Fforms_fun = Fforms_fun)
        }
        Times <- out$outcome.Times
        event <- out$outcome.event
        out <- out[names(out) != "outcome"]
        if (type == "ecdf") {
            r2 <- quantile(Times, probs = percentiles[2L], na.rm = TRUE)
            x_vals <- seq(0, r2, length.out = 500)
            rep_T <- apply(out$Times, 2L, function (x, x_vals) ecdf(x)(x_vals),
                           x_vals = x_vals)
            ss <- summary(survfit(Surv(Times, event) ~ 1), times = x_vals)
            F0 <- 1 - ss$surv
            F0_low <- 1 - ss$lower
            F0_upp <- 1 - ss$upper
            MISE <- mean(apply((rep_T - F0)^2, 2L, trapezoid_rule, x = x_vals))
            rootMISE <- round(sqrt(MISE), 5)
            if (is.null(xlab)) xlab <- "Event Times"
            if (is.null(ylab)) ylab <- "Empirical CDF"
            if (plot) {
                if (is.null(ylim)) {
                    yylim <- c(0, min(1, max(rep_T, F0_upp) + 0.1))
                } else yylim <- ylim
                matplot(x_vals, rep_T, type = "s", lty = lty_rep, lwd = lwd_rep,
                        col = col_rep, xlab = xlab, ylab = ylab, ylim = yylim, ...)
                title(main = main, line = line_main, cex.main = cex.main)
                lines(x_vals, F0, type = "s", lwd = lwd_obs, lty = lty_obs, col = col_obs)
                if (CI_ecdf != "none") {
                    lines(x_vals, F0_low, type = "s", lwd = lwd_obs, lty = 2, col = col_obs)
                    lines(x_vals, F0_upp, type = "s", lwd = lwd_obs, lty = 2, col = col_obs)
                }
                if (add_legend) {
                    if (!is.na(pos_legend[1L])) {
                        if (CI_ecdf != "none") {
                            legend(pos_legend[1L],
                                   legend = c("replicated data", "observed data", "95% CI"),
                                   lty = c(1, 1, 2), col = c(col_rep, col_obs, col_obs),
                                   bty = "n", cex = 0.9)
                        } else {
                            legend(pos_legend[1L],
                                   legend = c("replicated data", "observed data"),
                                   lty = c(1, 1), col = c(col_rep, col_obs),
                                   bty = "n", cex = 0.9)
                        }
                    }
                    if (!is.na(pos_legend[2L])) {
                        legend(pos_legend[2L], bty = "n",
                               legend = bquote(sqrt(MISE) == .(rootMISE)))
                    }
                }
            } else {
                plot_values <-
                    list(x_vals = x_vals, obs = F0, obs_low = F0_low,
                         obs_upp = F0_upp, rep = rep_T,
                         rootMISE = if (add_legend) rootMISE)
            }
        } else {
            DD <- data.frame(event = event)
            DD$S <- mapply(function (sims, times) ecdf(sims)(times),
                           sims = split(out$Times, row(out$Times)),
                           times = Times)
            KM <- survfit(Surv(S, event) ~ 1, data = DD)
            plot(KM, xlab = "empirical CDF at event times",
                 ylab = "Pr(U <= u)", main = main, fun = "event",
                 cex.main = cex.main, ...)
            xx <- seq(0, 1, length.out = 101)
            lines(xx, xx, lwd = 2, col = "red")
            if (add_legend) {
                legend("topleft", c("Kaplan-Meier transformed values",
                                    "CDF uniform"),
                       lty = 1, lwd = 2, col = c(1, 2), bty = "n")
            }
        }
    } else {
        list_of_jms <- inherits(object, 'list') && inherits(object[[1L]], 'jm') &&
            inherits(newdata, 'list')
        outY <- if (list_of_jms) {
            sims_per_fold <- if (is.null(params_mcmc)) {
                mapply2(simulate, object = object, newdata = newdata,
                        MoreArgs = list(process = "longitudinal",
                                        include_outcome = TRUE,
                                        random_effects = random_effects,
                                        return_random_effects = TRUE,
                                        seed = seed, nsim = nsim,
                                        Fforms_fun = Fforms_fun))
            } else {
                mapply2(simulate, object = object, newdata = newdata,
                        params_mcmc = params_mcmc,
                        MoreArgs = list(process = "longitudinal",
                                        include_outcome = TRUE,
                                        random_effects = random_effects,
                                        return_random_effects = TRUE,
                                        seed = seed, nsim = nsim,
                                        Fforms_fun = Fforms_fun))
            }
            bind(sims_per_fold)
        } else {
            simulate(object, nsim = nsim, newdata = newdata,
                     process = "longitudinal", include_outcome = TRUE,
                     random_effects = random_effects, params_mcmc = params_mcmc,
                     return_random_effects = TRUE,
                     seed = seed, Fforms_fun = Fforms_fun)
        }
        outT <- if (inherits(object, 'list') && inherits(object[[1L]], 'jm') &&
                   inherits(newdata, 'list')) {
            sims_per_fold <- if (is.null(params_mcmc)) {
                mapply2(simulate, object = object, newdata = newdata,
                        MoreArgs = list(process = "event",
                                        include_outcome = TRUE,
                                        Fforms_fun = Fforms_fun,
                                        random_effects = attr(outY[[1L]], "RE"),
                                        seed = seed, nsim = nsim))
            } else {
                mapply2(simulate, object = object, newdata = newdata,
                        params_mcmc = params_mcmc,
                        MoreArgs = list(process = "event",
                                        include_outcome = TRUE,
                                        Fforms_fun = Fforms_fun,
                                        random_effects = attr(outY[[1L]], "RE"),
                                        seed = seed, nsim = nsim))
            }
            bind(sims_per_fold)
        } else {
            simulate(object, nsim = nsim, newdata = newdata,
                     process = "event", seed = seed, include_outcome = TRUE,
                     random_effects = attr(outY[[1L]], "RE"),
                     params_mcmc = params_mcmc, Fforms_fun = Fforms_fun)
        }
        Y <- outY$outcome
        n_outcomes <- length(Y)
        index <- seq_len(n_outcomes)
        if (outcomes < Inf) index <- index[index %in% outcomes]
        nindex <- length(index)
        id <- lapply(Y, function (x) attr(x, "id"))
        times <- lapply(Y, function (x) attr(x, "times"))
        outY <- outY[names(outY) != "outcome"]
        Time <- outT$outcome.Times
        event <- outT$outcome.event
        outT <- outT[names(outT) != "outcome"]
        for (j in seq_len(n_outcomes)) {
            id_j <- id[[j]]
            times_j <- times[[j]]
            f <- function (t, T) {t[t > T] <- NA; t}
            for (m in seq_len(ncol(outT$Times))) {
                nas <- is.na(unlist(mapply2(f, t = split(times_j, id_j),
                                            T = outT$Times[, m])))
                outY[[j]][nas, m] <- NA_real_
            }
        }
        association <- function (Time, event, Y, times, id, Uno = FALSE) {
            n_outcomes <- length(Y)
            YY <- mapply2(split, Y, id)
            ttimes <- mapply2(split, times, id)
            unq_eventTimes <- c(0, sort(unique(Time[event == 1])))
            unq_eventTimes <- unq_eventTimes[unq_eventTimes < quantile(unq_eventTimes, 0.95)]
            f <- function (y, t, t0) {
                if (all(t > t0)) NA else y[max(which((t - t0) <= 0), na.rm = TRUE)]
            }
            Cs <- matrix(0, length(unq_eventTimes), n_outcomes)
            for (i in seq_along(unq_eventTimes)) {
                t0 <- unq_eventTimes[i]
                for (j in seq_len(n_outcomes)) {
                    DF <- data.frame(Time = Time, event = event)
                    DF[["Y"]] <-
                        mapply(f, y = YY[[j]], t = ttimes[[j]],
                               MoreArgs = list(t0 = t0))
                    DF <- DF[!is.na(DF[["Y"]]), ]
                    DF$lp <- coxph(Surv(Time, event) ~ Y, data = DF)$linear.predictors
                    Cs[i, j] <- if (Uno) {
                        survC1::Est.Cval(cbind(DF$Time, DF$event, DF$lp),
                                         max(Time) - 0.001, nofit = TRUE)$Dhat
                    } else {
                        concordance(Surv(Time, event) ~ lp, data = DF,
                                    reverse = TRUE)$concordance
                    }
                }
            }
            lapply(seq_len(n_outcomes), function (j)
                cbind(unq_eventTimes, C = Cs[, j]))
        }
        loess.smooth2 <- function (x, y) {
            loess.smooth(x, y, degree = 1, span = 0.75,
                         family = "gaussian", evaluation = 200)
        }
        C <- association(Time, event, Y, times, id, TRUE)
        Obs <- mapply2(loess.smooth2, x = lapply(C, function (c) c[, 1L]),
                       y = lapply(C, function (c) c[, 2L]))
        if (CI_loess) {
            low <- upp <- vector("list", n_outcomes)
            for (j in seq_len(n_outcomes)) {
                loess_fit <-
                    loess(y ~ x, data = data.frame(x = C[[j]][, 1L], y = C[[j]][, 2L]),
                          degree = 1L, span = 0.75, family = "gaussian")
                preds <- predict(loess_fit, data.frame(x = Obs[[j]]$x), se = TRUE)
                low[[j]] <- preds$fit + qt(0.025, preds$df) * preds$se.fit
                upp[[j]] <- preds$fit + qt(0.975, preds$df) * preds$se.fit
            }
        }
        Rep <- vector("list", ncol(outY[[1L]]))
        for (m in seq_along(Rep)) {
            C_m <- association(outT$Times[, m], outT$event[, m],
                               lapply(outY, function (y) y[, m]),
                               times, id)
            Rep[[m]] <-
                mapply2(loess.smooth2, x = lapply(C_m, function (c) c[, 1L]),
                        y = lapply(C_m, function (c) c[, 2L]))
        }
        MISE <- numeric(n_outcomes)
        for (j in seq_len(n_outcomes)) {
            Obs_j <- Obs[[j]]
            Rep_j <- lapply(Rep, function (x) x[[j]])
            MISE[j] <- mean(sapply(Rep_j, function (r, o)
                trapezoid_rule((r$y - o$y)^2, o$x), o = Obs_j))
        }
        rootMISE <- round(sqrt(MISE), 5)
        if (plot) {
            if (is.null(xlab)) xlab <- "Event Times"
            if (is.null(ylab)) ylab <- "Concordance"
            xlab <- rep(xlab, length.out = nindex)
            ylab <- rep(ylab, length.out = nindex)
            main <- rep(main, length.out = nindex)
            for (j in index) {
                jj <- match(j, index)
                rx <- range(c(sapply(Rep, function (loe) loe[[j]]$x), Obs[[j]]$x))
                ry <- range(c(sapply(Rep, function (loe) loe[[j]]$y), Obs[[j]]$y))
                if (CI_loess) {
                    ry <- range(ry, low[[j]], upp[[j]])
                }
                if (is.null(ylim)) {
                    yylim <- ry
                } else yylim <- ylim
                plot(rx, ry, type = "n", xlab = xlab[jj], ylab = ylab[jj],
                     main = main[jj], cex.main = cex.main, ylim = yylim, ...)
                for (i in seq_along(Rep)) {
                    lines(Rep[[i]][[j]], col = col_rep, lty = lty_rep, lwd = lwd_rep)
                }
                lines(Obs[[j]], col = col_obs, lty = lty_obs, lwd = lwd_obs)
                if (CI_loess) {
                    lines(Obs[[j]]$x, low[[j]], lwd = lwd_obs, lty = 2, col = col_obs)
                    lines(Obs[[j]]$x, upp[[j]], lwd = lwd_obs, lty = 2, col = col_obs)
                }
                if (add_legend) {
                    if (!is.na(pos_legend[1L])) {
                        if (CI_loess) {
                            legend(pos_legend[1L],
                                   legend = c("replicated data", "observed data", "95% CI"),
                                   lty = c(1, 1, 2), col = c(col_rep, col_obs, col_obs),
                                   bty = "n", cex = 0.9)
                        } else {
                            legend(pos_legend[1L],
                                   legend = c("replicated data", "observed data"),
                                   lty = c(1, 1), col = c(col_rep, col_obs),
                                   bty = "n", cex = 0.9)
                        }
                    }
                    if (!is.na(pos_legend[2L])) {
                        legend(pos_legend[2L], bty = "n",
                               legend = bquote(sqrt(MISE) == .(rootMISE[j])))
                    }
                }
            }
        } else {
            plot_values <-
                list(x_vals = Obs, obs = Obs,
                     obs_low = if (CI_loess) low,
                     obs_upp = if (CI_loess) upp, rep = Rep,
                     rootMISE = if (add_legend) NA)
        }
    }
    if (!plot) return(plot_values)
}

# nlme::lme() / survival::coxph() do not dispatch on the class of the `data` argument.
# We define our own S3 generics that dispatch on `data` (2nd argument),
# so that lme.sliced_data() / mixed_model.sliced_data() / coxph.sliced_data()
# are used when `data` is `sliced_data`, while nlme::lme() /
# GLMMadaptive::mixed_model() / survival::coxph() are used for `data.frame`
lme <- function (fixed, data, ...) UseMethod ("lme", data)
mixed_model <- function (fixed, data, ...) UseMethod ("mixed_model", data)
coxph <- function (formula, data, ...) UseMethod ("coxph", data)
jm <- function (Surv_object, Mixed_objects, time_var, recurrent = FALSE,
                functional_forms = NULL, which_independent = NULL,
                base_hazard = NULL, data_Surv = NULL, id_var = NULL,
                priors = NULL, control = NULL, ...) {
  UseMethod ("jm", Surv_object)
}

# We cannot implement coxph.data.frame() as a simple forwarder like
# coxph.data.frame <- function (formula, data, ...) {
#   survival::coxph(formula = formula, data = data, ...)
# }
# because the resulting coxph object will typically store the call with
# `object$call$data` as `data`, where `data` is just the argument name of this
# method, not the original expression supplied by the user (e.g., `pbc2.id`).
# Later when jm() tries to recover the survival dataset via
# `eval(Surv_object$call$data, parent.frame())`, if `object$call$data` is `data`,
# jm() may end up evaluating whatever object is named `data` in the environment.
# For this reason, we reconstruct the call using match.call() and evaluate it in the
# parent frame, so object$call$data retains the original expression (e.g., `pbc2.id`).
coxph.default <- function (formula, data, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- quote(survival::coxph)
  eval(mc, parent.frame())
}

lme.default <- function (fixed, data, random, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- quote(nlme::lme)
  eval(mc, parent.frame())
}

mixed_model.default <- function (fixed, data, random, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- quote(GLMMadaptive::mixed_model)
  eval(mc, parent.frame())
}

slapply <- function (X, FUN, ..., parallel = TRUE, ncores = NULL,
                     pkgs = NULL, backend = c("auto", "multicore", "snow")) {
  backend <- match.arg(backend)

  if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)

  if (!parallel || ncores <= 1) return(lapply(X, FUN, ...))

  if (backend == "auto") {
    # multicore: parallel::mclapply() uses OS forking (Unix/macOS/Linux only).
    #            Fast/low overhead because workers are forked copies of the
    #            current R process. Downside: forking can be fragile if the
    #            current R session has loaded "native" compiled code (C/C++/Fortran)
    #            that uses threads internally (e.g., %*%, chol()). In those cases,
    #            forked workers may hang/crash or behave unpredictably.
    # snow: PSOCK socket cluster, parallel::makeCluster(), starts fresh R worker
    #       sessions. Slower startup and more data transfer (packages must be
    #       loaded on workers), but generally the most robust and cross-platform
    #       option (Windows/Linux/macOS).
    backend <- if (.Platform$OS.type == "windows") "snow" else "multicore"
  }

  if (backend == "multicore") {
    return(parallel::mclapply(X, FUN, ..., mc.cores = ncores))
  }

  # backend == "snow"
  cl <- parallel::makeCluster(ncores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  if (!is.null(pkgs)) { # load packages on workers
    pkgs <- as.character(pkgs)
    load_pkgs <- function (pkgs) {
      for (p in pkgs) {
        suppressPackageStartupMessages(library (p, character.only = TRUE))
      }
      NULL
    }
    environment(load_pkgs) <- baseenv()
    parallel::clusterCall(cl, load_pkgs, pkgs = pkgs)
  }

  parallel::parLapply(cl, X, FUN, ...)
}

coxph.sliced_data <- function (formula, data, ...,
                               parallel_out = TRUE, cores = NULL) {
  dots <- list(...)
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)

  fit  <- function (dt, formula, dots) {
    args <- c(list(formula = formula, data = dt), dots)
    do.call(survival::coxph, args)
  }

  environment(fit) <- new.env(parent = baseenv()) # Detach fit from the parent call frame so it doesn't accidentally serialize large objects (i.e., the full `data` list) when sent to "snow" workers. print(ls(environment(fit)))

  fits <- slapply(X = data, FUN = fit,
                  formula = formula, dots = dots,
                  parallel = parallel_out, ncores = cores,
                  pkgs = "survival")

  class(fits) <- c("sliced_coxph", class(fits))
  fits
}

lme.sliced_data <- function (fixed, data, random, ...,
                             parallel_out = TRUE, cores = NULL) {
  dots <- list(...)
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)

  fit <- function (dt, fixed, random, dots) {
    args <- c(list(fixed = fixed, data = dt, random = random), dots)
    do.call(nlme::lme, args)
  }

  environment(fit) <- new.env(parent = baseenv()) # Detach fit from the parent call frame so it doesn't accidentally serialize large objects (i.e., the full `data` list) when sent to "snow" workers. print(ls(environment(fit)))

  fits <- slapply(X = data, FUN = fit,
                  fixed = fixed, random = random, dots = dots,
                  parallel = parallel_out, ncores = cores,
                  pkgs = "nlme")

  class(fits) <- c("sliced_lme", class(fits))
  fits
}

mixed_model.sliced_data <- function (fixed, data, random, ...,
                                     parallel_out = TRUE, cores = NULL) {
  dots <- list(...)
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)

  bool1 <- parallel_out && cores > 1L
  bool2 <- identical(dots$optimizer, "optimParallel") ||
    (is.list(dots$control) && identical(dots$control$optimizer, "optimParallel")) # identical(NULL, "optimParallel") -> FALSE
  if (bool1 && bool2) {
    warning(
      "'optimParallel' requested, but mixed_model.sliced_data() is already parallel across subsamples; ",
      "switching optimizer to 'optim' to avoid nested parallelism. ",
      "To use 'optimParallel', run mixed_model.sliced_data(parallel_out = FALSE).",
      call. = FALSE
    )
    if (identical(dots$optimizer, "optimParallel")) dots$optimizer <- "optim"
    if (is.list(dots$control) && identical(dots$control$optimizer, "optimParallel")) {
      dots$control$optimizer <- "optim"
    }
  }

  fit <- function (dt, fixed, random, dots) {
    args <- c(list(fixed = fixed, random = random, data = dt), dots)
    do.call(GLMMadaptive::mixed_model, args)
  }

  environment(fit) <- new.env(parent = baseenv()) # Detach fit from the parent call frame so it doesn't accidentally serialize large objects (i.e., the full `data` list) when sent to "snow" workers. print(ls(environment(fit)))

  fits <- slapply(X = data, FUN = fit,
                  fixed = fixed, random = random, dots = dots,
                  parallel = parallel_out, ncores = cores,
                  pkgs = "GLMMadaptive", backend = "snow")

  class(fits) <- c("sliced_MixMod", class(fits))
  fits
}

jm.sliced_coxph <- function (Surv_object, Mixed_objects, time_var, ...,
                             parallel_out = TRUE, cores = NULL) {
  n_slices <- length(Surv_object)
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)

  dots <- list(...)
  ctrl <- dots$control
  if (is.null(ctrl)) ctrl <- list()

  # n_chains: control$n_chains > n_chains > default
  if (!is.null(ctrl$n_chains)) {
    n_chains <- ctrl$n_chains
  } else if (!is.null(dots$n_chains)) {
    n_chains <- dots$n_chains
  } else {
    n_chains <- 3L
  }

  # inner jm() backend (chains)
  # parallel: control$parallel > parallel > auto
  if (is.null(ctrl$parallel) && !is.null(dots$parallel)) ctrl$parallel <- dots$parallel
  if (is.null(ctrl$parallel)) {
    ctrl$parallel <- if (.Platform$OS.type == "windows") "snow" else "multicore"
  }
  if (identical(ctrl$parallel, "multicore") && .Platform$OS.type == "windows") {
    warning("control$parallel='multicore' is not available on Windows; using 'snow'.",
            call. = FALSE)
    ctrl$parallel <- "snow"
  }

  inner_ncores <- max(1L, min(n_chains, cores))
  ctrl$cores <- inner_ncores
  ctrl$n_chains <- n_chains

  dots[c("n_chains", "cores", "parallel")] <- NULL # Avoid passing control-like args twice.
  dots$control <- ctrl

  # outer cores (across slices)
  outer_ncores <- min(n_slices, max(1L, floor(cores / inner_ncores)))

  # Optional: on Windows, avoid nested PSOCK clusters by default
  # if (.Platform$OS.type == "windows" && inner_ncores > 1L && identical(ctrl$parallel, "snow")) {
  #   outer_ncores <- 1L
  # }

  # Detect whether Mixed_objects is:
  # i)  a sliced object: list of lme/MixMod fits (one per slice)
  # ii) a list of sliced objects: list(outcome1_sliced, outcome2_sliced, ...) -> transpose
  if (!inherits(Mixed_objects[[1]], c("lme", "MixMod"))) {
    Mixed_objects <- lapply(seq_len(n_slices),
                            function (i) lapply(Mixed_objects, `[[`, i))
  }

  jobs <- Map(function(S, M) list(Surv_object = S, Mixed_objects = M), # Each worker receives only slice-specific inputs.
              Surv_object, Mixed_objects)

  fit <- function (job, time_var, dots) {
    args <- c(list(Surv_object = job$Surv_object,
                   Mixed_objects = job$Mixed_objects,
                   time_var = time_var),
              dots)
    do.call(JMbayes2::jm, args)
  }

  environment(fit) <- new.env(parent = baseenv())

  fits <- slapply(X = jobs, FUN = fit,
                  time_var = time_var, dots = dots,
                  parallel = parallel_out, ncores = outer_ncores,
                  pkgs = "JMbayes2", backend = "snow")

  class(fits) <- c("sliced_jm", class(fits))
  fits
}

slicer <- function (n_slices, id_var, data_long, data_surv, seed = 123L) {
  stopifnot(is.numeric(n_slices), length(n_slices) == 1, n_slices >= 1,
            is.character(id_var))
  if (!id_var %in% names(data_long)) {
    stop(paste0("'", id_var, "' not found in data_long."))
  }
  if (!id_var %in% names(data_surv)) {
    stop(paste0("'", id_var, "' not found in data_surv."))
  }
  ids_unq <- unique(c(as.character(data_long[[id_var]]),
                      as.character(data_surv[[id_var]])))
  if (!is.null(seed)) set.seed(seed)
  grp <- ((seq_along(ids_unq) - 1) %% n_slices) + 1 # assign each ID a group number 1...n_slices in round-robin order (1,2,...,n_slices,1,2,...)
  ids_slc <- split(sample(ids_unq), grp)
  long <- lapply(ids_slc, function (ids) data_long[data_long[[id_var]] %in% ids, ])
  surv <- lapply(ids_slc, function (ids) data_surv[data_surv[[id_var]] %in% ids, ])
  class(long) <- c("sliced_data", class(long))
  class(surv) <- c("sliced_data", class(surv))
  list(long = long, surv = surv)
}

consensus <- function (object, parm,
                       method = c("union", "equal_weight", "var_weight"),
                       seed = 123L) {
  stopifnot(inherits(object, "sliced_jm"), is.character(parm), length(parm) >= 1)
  method <- match.arg(method)
  get_cons <- function (fits, parm, method) {
    slice_mats <- lapply(fits, function (x) do.call(rbind, x$mcmc[[parm]])) # list of [iter, p] matrices (one per slice)
    if (method == "union") {
      return (list(draws = do.call(rbind, slice_mats), # [iter * slice, p]
                   weights = NULL))
    }
    draws <- simplify2array(slice_mats) # [iter, p, slice]
    draws <- apply(draws, c(2, 3), sample)
    ns <- dim(draws)[3]
    snames <- paste0("S", seq_len(ns)) # slice names
    pnames <- dimnames(draws)[[2]]
    if (method == "equal_weight") {
      w <- matrix(1 / ns, nrow = length(pnames), ncol = ns,
                  dimnames = list(pnames, snames))
      return (list(draws = apply(draws, c(1, 2), mean), # [iter, p]
                   weights = w))
    }
    # var_weight
    vars <- apply(draws, c(2, 3), var) # [p, slice]
    w <- 1 / pmax(vars, .Machine$double.eps) # [p, slice]
    w <- sweep(w, 1, rowSums(w), "/") # normalize
    dimnames(w) <- list(pnames, snames)
    w_draws <- sweep(draws, c(2, 3), w, "*") # weighted[iter, p, slice] <- draws[iter, p, slice] * w[p, slice]
    list(draws = apply(w_draws, c(1, 2), sum), weights = w)
  }
  if (!is.null(seed)) set.seed(seed)
  cons_out <- lapply(parm, function (p) get_cons(object, p, method))
  names(cons_out) <- parm
  cons_draws <- lapply(cons_out, `[[`, "draws")
  cons_wts   <- lapply(cons_out, `[[`, "weights")
  n_draws <- vapply(cons_draws, nrow, integer(1))
  summarise_draws <- function (draws_mat) {
    cbind(
      Mean   = colMeans(draws_mat),
      StDev  = apply(draws_mat, 2, sd),
      `2.5%`  = apply(draws_mat, 2, quantile, probs = 0.025, names = FALSE),
      `97.5%` = apply(draws_mat, 2, quantile, probs = 0.975, names = FALSE),
      P      = apply(draws_mat, 2, function (x) 2 * min(mean(x > 0), mean(x < 0)))
    )
  }
  cons_sum <- lapply(cons_draws, summarise_draws)
  res <- list(method = method, parm = parm, n_splits = length(object),
              draws = cons_draws, weights = cons_wts, n_draws = n_draws,
              summary = cons_sum, seed = seed)
  class(res) <- "consensus_jm"
  res
}

print.consensus_jm <- function (x, digits = 4, ...) {
  cat("\nConsensus summary (", x$n_splits, " splits)\n", sep = "")
  if (identical(x$method, "union")) {
    cat("Method (union): concatenated draws across splits (no averaging).\n")
  } else if (identical(x$method, "equal_weight")) {
    cat("Method (equal_weight): iteration-wise simple average across splits.\n")
  } else if (identical(x$method, "var_weight")) {
    cat("Method (var_weight): iteration-wise weighted average across splits using inverse-variance weights.\n")
  }
  for (nm in names(x$summary)) {
    cat("\nParameter block: ", nm,
        " (", x$n_draws[[nm]], " draws)\n", sep = "")
    tab <- x$summary[[nm]]
    w <- x$weights[[nm]]
    if (!is.null(w)) {
      w <- w[rownames(tab), , drop = FALSE]
      colnames(w) <- paste0("w_", colnames(w))
      tab <- cbind(tab, w)
    }
    print(round(tab, digits))
  }
  invisible(x)
}

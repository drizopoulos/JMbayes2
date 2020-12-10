traceplot <- function (object, ...) UseMethod("traceplot")

#' Trace plot of MCMC output for Joint Models
#' 
#' Plots the evolution of the estimated parameter vs. iterations in a 
#' fitted joint model.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param parm  either \code{"all"} or one specific joint model parameter of interest.
#' @param ... further arguments passed to \code{\link[coda]{traceplot}}.
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link[coda]{traceplot}}, \code{\link{ggtraceplot}}, 
#' \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fit
#' fit_lme <- lme(sqrt(CD4) ~ obstime * drug, random = ~ 1 + obstime | patient, 
#'                data = aids)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(Time, death) ~ drug, data = aids.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, fit_lme, time_var = "obstime")
#' 
#' # trace plot for the fixed effects in the linear mixed submodel 
#' traceplot(fit_jm, parm = "betas")
#' }
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

#' Gelman and Rubin's Convergence Diagnostic for Joint Models
#' 
#' Calculates the potential scale reduction factor for the estimated parameters
#' in a fitted joint model, together with the upper confidence limits.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param parm  either \code{"all"} or one specific joint model parameter of interest.
#' @param ... further arguments passed to \code{\link[coda]{gelman.diag}}.
#' @return A list of \code{gelman.diag} objects. An object of class 
#' \code{gelman.diag} is a list with the elements:
#' \tabular{ll}{
#' \code{psrf} \tab a list containing the point estimates of the potential 
#' scale reduction factor (labelled \code{Point est.}) and their upper 
#' confidence limits (labelled \code{Upper C.I.}). \cr
#' \code{mpsrf} \tab the point estimate of the multivariate potential scale 
#' reduction factor. This is NULL if the parameter is univariate.
#' }
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link[coda]{gelman.diag}}, \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fit
#' fit_lme <- lme(sqrt(CD4) ~ obstime * drug, 
#' random = ~ 1 + obstime | patient, data = aids)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(Time, death) ~ drug, data = aids.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, fit_lme, time_var = "obstime")
#' 
#' # joint model convergence diagnostic
#' gelman_diag(fit_jm, "all")
#' }
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

densplot <- function (object, ...) UseMethod("densityplot")

#' Probability Density Plot for Joint Models
#' 
#' Plots the density estimate for the estimated parameters
#' in a fitted joint model.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param parm  either \code{"all"} or one specific joint model parameter of interest.
#' @param ... further arguments passed to \code{\link[coda]{densplot}}.
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link[coda]{densplot}}, \code{\link{ggdensityplot}}, 
#' \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fit
#' fit_lme <- lme(sqrt(CD4) ~ obstime * drug, random = ~ 1 + obstime | patient,
#'                data = aids)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(Time, death) ~ drug, data = aids.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, fit_lme, time_var = "obstime")
#' 
#' # density plot for the fixed effects in all linear mixed submodels 
#' densplot.jm(fit_jm, parm = "betas")
#' }
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

#' Cumulative Quantile Plot for Joint Models
#' 
#' Plots the evolution of the sample quantiles vs. iterations in a fitted joint 
#' model.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param parm  either \code{"all"} or one specific joint model parameter of interest.
#' @param ... further arguments passed to \code{\link[coda]{cumuplot}}.
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link[coda]{cumuplot}}, \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fit
#' fit_lme <- lme(sqrt(CD4) ~ obstime * drug, random = ~ 1 + obstime | patient, 
#'                data = aids)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(Time, death) ~ drug, data = aids.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, fit_lme, time_var = "obstime")
#' 
#' # cumulative quantile plot for the fixed effects in all linear mixed submodels 
#' cumuplot(fit_jm, parm = "betas")
#' }
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

#' Fixed Effects Estimates for Survival Submodel in Joint Models
#' 
#' Extracts estimated fixed effects for the event process from a fitted joint model.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param ... additional arguments; currently none is used.
#' @return A list with the elements: 
#' \tabular{ll}{
#' \code{gammas} \tab estimated baseline fixed effects. \cr
#' \code{association} \tab estimated association parameters.   
#' }
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link{fixef}}, \code{\link{ranef}}, \code{\link{jm}}.
#' @examples 
#' \dontrun{
#' # linear mixed model fit
#' fit_lme <- lme(sqrt(CD4) ~ obstime * drug, random = ~ 1 + obstime | patient, 
#'                data = aids)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(Time, death) ~ drug, data = aids.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, fit_lme, time_var = "obstime")
#' 
#' # fixed effects for the event process 
#' coef(fit_jm)
#' }
coef.jm <- function (object, ...) {
    gammas <- object$statistics$Mean[["gammas"]]
    if (is.null(gammas)) object$statistics$Mean[["alphas"]] else
        list("gammas" = gammas,
             "association" = object$statistics$Mean[["alphas"]])
}

#' Fixed Effects Estimates for Linear Mixed Submodels in Joint Models
#' 
#' Extracts estimated fixed effects for the longitudinal processes from a fitted joint model.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param outcome the linear mixed submodel to extract the estimated fixed effects. 
#' If greater than the total number of linear mixed submodels, extracts from all 
#' of them.
#' @param ... additional arguments; currently none is used.
#' @return A numeric vector of the estimated fixed effects for the 
#' \code{outcome} selected. If \code{outcome} is greater than the number of 
#' linear mixed submodels, returns a list of numeric vectors for all outcomes.
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link{coef}}, \code{\link{ranef}}, \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fit
#' fit_lme1 <- lme\(log(serBilir) ~ year:sex + age,
#'                 random = ~ year | id, data = pbc2)
#' 
#' fit_lme2 <- lme(prothrombin ~ sex, 
#'                 random = ~ year | id, data = pbc2)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, list(fit_lme1, fit_lme2), time_var = "year")
#' 
#' # fixed effects for all linear mixed submodels 
#' fixef(fit_jm, outcome = 3)
#' }
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

#' Random Effects Estimates for Linear Mixed Submodels in Joint Models
#' 
#' Extracts estimated random effects for the longitudinal processes from a fitted joint model.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param outcome the linear mixed submodel to extract the estimated fixed effects.
#' If greater than the total number of linear mixed submodels, extracts from all 
#' of them.
#' @param post_vars logical; if TRUE returns the variance of the posterior distribution.
#' @param ... additional arguments; currently none is used.
#' @return A numeric matrix with rows denoting the individuals and columns the random 
#' effects. If \code{postVar = TRUE}, the numeric matrix has an extra attribute 
#' “postVar".
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link{coef}}, \code{\link{fixef}}, \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fits
#' fit_lme1 <- lme(log(serBilir) ~ year:sex + age, random = ~ year | id, 
#'                 data = pbc2)
#' 
#' fit_lme2 <- lme(prothrombin ~ sex, random = ~ year | id, data = pbc2)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, list(fit_lme1, fit_lme2), time_var = "year")
#' 
#' # random effects from all linear mixed submodels 
#' ranef(fit_jm)
#' }
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

#' Joint Model Terms
#' 
#' Extracts the terms objects from a fitted joint model.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param process which submodel to extract the terms, i.e.,
#' \tabular{ll}{
#' \code{longitudinal} \tab linear mixed model(s). \cr
#' \code{survival} \tab survival model.   
#' }
#' @param type which effects, fixed or random, to select, when 
#' \code{process = "longitudinal"}.
#' @param ... additional arguments; currently none is used.
#' @return If \code{process = "longitudinal"}, a list of the terms object(s) for
#'  the linear mixed model(s). If \code{process = "event"}, the terms object 
#'  for the survival model.
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link{model.frame}}, \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fit
#' fit_lme <- lme(log(serBilir) ~ year * sex, random = ~ year | id, data = pbc2)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, fit_lme, time_var = "year")
#' 
#' # fixed effects terms in the linear mixed model
#' terms(fit_jm, process = "longitudinal", type = "random")
#' }
terms.jm <- function (object, process = c("longitudinal", "event"),
                      type = c("fixed", "random"), ...) {
    process <- match.arg(process)
    type <- match.arg(type)
    combo <- paste(process, type, sep = "_")
    switch(combo,
           "longitudinal_fixed" = object$model_info$terms$terms_FE,
           "longitudinal_random" = object$model_info$terms$terms_RE,
           "event_fixed" = , "event_random" = object$model_info$terms$terms_Surv)
}

#' Model.frame method for Joint Models
#' 
#' Creates the model frame from a fitted joint model.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param process which submodel to recreate the model frame, i.e.,
#' \tabular{ll}{
#' \code{longitudinal} \tab linear mixed model(s). \cr
#' \code{survival} \tab survival model.   
#' }
#' @param type which effects, fixed or random, to select, when 
#' \code{process = longitudinal}.
#' @param ... additional arguments; currently none is used.
#' @return If \code{process = "longitudinal"}, a list of the model frames used 
#' in the linear mixed model(s). If \code{process = "event"}, the model frame 
#' used in the survival model.
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link[stats]{model.frame}}, \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fit
#' fit_lme1 <- lme(log(serBilir) ~ year:sex + age, random = ~ year | id, 
#'                 data = pbc2)
#'
#' fit_lme2 <- lme(prothrombin ~ sex, random = ~ year | id, data = pbc2)
#'
#' # cox model fit
#' fit_cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)
#'
#' # joint model fit
#' fit_jm <- jm(fit_cox, list(fit_lme1, fit_lme2), time_var = "year")
#'
#' # model frame for the fixed effects terms in the linear mixed submodels
#' model.frame(fit_jm, process = "longitudinal", type = "fixed")
#' }
model.frame.jm <- function (object, process = c("longitudinal", "event"),
                            type = c("fixed", "random"), ...) {
    process <- match.arg(process)
    type <- match.arg(type)
    combo <- paste(process, type, sep = "_")
    switch(combo,
           "longitudinal_fixed" = object$model_info$frames$mf_FE,
           "longitudinal_random" = object$model_info$frames$mf_RE,
           "event_fixed" = , "event_random" = object$model_info$frames$mf_Surv)
}

#' Design Matrices for Linear Mixed Submodels in Joint Models
#' 
#' Creates the design matrices for linear mixed submodels from a fitted joint 
#' model.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param ... additional arguments; currently none is used.
#' @return A list of the design matrices for the linear mixed submodels.
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link[stats]{model.matrix}}, \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fit
#' fit_lme1 <- lme(log(serBilir) ~ year:sex + age, random = ~ year | id, 
#'                 data = pbc2)
#'
#' fit_lme2 <- lme(prothrombin ~ sex, random = ~ year | id, 
#'                 data = pbc2)
#'
#' # cox model fit
#' fit_cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)
#'
#' # joint model fit
#' fit_jm <- jm(fit_cox, list(fit_lme1, fit_lme2), time_var = "year")
#'
#' # linear mixed models design matrices
#' model.matrix(fit_jm)
#' }
model.matrix.jm <- function (object, ...) {
    tr <- terms(object)
    mf <- model.frame(object)
    if (is.data.frame(mf)) {
        model.matrix(tr, mf)
    } else {
        mapply(model.matrix.default, object = tr, data = mf, SIMPLIFY = FALSE)
    }
}

#' Family Objects for Joint Models
#' 
#' Extracts the error distribution and link function used in the linear
#' mixed submodel(s) from a fitted joint model.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param ... additional arguments; currently none is used.
#' @return A list of \code{family} objects.
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link[stats]{family}}, \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fits
#' fit_lme <- lme(log(serBilir) ~ year * sex, random = ~ year | id, 
#'                data = pbc2)
#'
#' fit_mm <- mixed_model(ascites ~ year * sex, random = ~ 1 | id, 
#'                       family = binomial(), data = pbc2)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, list(fit_lme, fit_mm), time_var = "year")
#' 
#' # family objects from all linear mixed submodels
#' family(fit_jm)
#' }
family.jm <- function (object, ...) {
    object$model_info$families
}

# ggplot mcmc diagnostics need ggplot2 and gridExtra
ggtraceplot <- function (object, ...) UseMethod("ggtraceplot")


#' Trace plot of MCMC output for Joint Models using ggplot2
#' 
#' Plots the evolution of the estimated parameter vs. iterations in a 
#' fitted joint model using ggplot2.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param parm  a character string with either the value \code{"all"} or one specific joint model group of parameters of interest. Defaults to \code{'all'}. \cr
#' @param size  the width of the traceplot line in mm. Defaults to 1.
#' @param alpha the opacity level of the traceplot line. Defaults to 0.8.
#' @param theme a character string specifying the color theme to be used. Possible options are \code{'standard'}, \code{'catalog'}, \code{'metro'}, \code{'pastel'}, \code{'beach'}, \code{'moonlight'}, \code{'goo'}, \code{'sunset'}. \cr
#' @param grid  logical; defaults to \code{FALSE}. If \code{TRUE} the plots are returned in grids split over multiple pages. For more details see the documentation for \code{\link[gridExtra:arrangeGrob]{gridExtra::marrangeGrob()}}. \cr  
#' @param gridrows number of rows per page for the grid. Only relevant when using \code{grid = TRUE}. Defaults to 3.
#' @param gridcols number of columns per page for the grid. Only relevant when using \code{grid = TRUE}. Defaults to 1. 
#' @param ... additional arguments; currently none is used.
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link{traceplot}}, 
#' \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fit
#' fit_lme <- lme(sqrt(CD4) ~ obstime * drug, random = ~ 1 + obstime | patient, 
#'                data = aids)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(Time, death) ~ drug, data = aids.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, fit_lme, time_var = "obstime")
#' 
#' # trace plot for the fixed effects in the linear mixed submodel 
#' ggtraceplot(fit_jm, parm = "betas")
#' ggtraceplot(fit_jm, parm = "betas", grid = TRUE)
#' }

# traceplot with ggplot
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
                geom_line(aes(x = iteration, y = value, color = chain), size = size, alpha = alpha) +
                ggtitle(paste('Traceplot of ', unique(ggdata$parm)[i])) +
                theme_bw() + theme(plot.title = element_text(hjust=0.5)) +
                scale_color_manual(values = ggcolthemes[[coltheme]]) +
                guides(color = guide_legend(override.aes = list(alpha = 1)))
        }
        marrangeGrob(grobs = gplots, nrow = gridrows, ncol = gridcols)
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

#' Probability Density Plot for Joint Models using ggplot2
#' 
#' Plots the density estimate for the estimated parameters
#' in a fitted joint model using ggplot2.
#' 
#' @param object an object inheriting from class \code{"jm"}.
#' @param parm  a character string with either the value \code{"all"} or one specific joint model group of parameters of interest. Defaults to \code{'all'}. \cr
#' @param size  the width of the density outline in mm. Defaults to 1.
#' @param alpha the opacity level of the density plot. Defaults to 0.6.
#' @param theme a character string specifying the color theme to be used. Possible options are \code{'standard'}, \code{'catalog'}, \code{'metro'}, \code{'pastel'}, \code{'beach'}, \code{'moonlight'}, \code{'goo'}, \code{'sunset'}. \cr
#' @param grid  logical; defaults to \code{FALSE}. If \code{TRUE} the plots are returned in grids split over multiple pages. For more details see the documentation for \code{\link[gridExtra:arrangeGrob]{gridExtra::marrangeGrob()}}. \cr  
#' @param gridrows number of rows per page for the grid. Only relevant when using \code{grid = TRUE}. Defaults to 3.
#' @param gridcols number of columns per page for the grid. Only relevant when using \code{grid = TRUE}. Defaults to 1. 
#' @param ... additional arguments; currently none is used.
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link{densplot}}, 
#' \code{\link{jm}}.
#' @examples
#' \dontrun{
#' # linear mixed model fit
#' fit_lme <- lme(sqrt(CD4) ~ obstime * drug, random = ~ 1 + obstime | patient, 
#'                data = aids)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(Time, death) ~ drug, data = aids.id)
#' 
#' # joint model fit
#' fit_jm <- jm(fit_cox, fit_lme, time_var = "obstime")
#' 
#' # trace plot for the fixed effects in the linear mixed submodel 
#' ggdensityplot(fit_jm, parm = "betas")
#' ggdensityplot(fit_jm, parm = "betas", grid = TRUE)
#' }

# density plot with ggplot
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
                geom_density(aes(x = value, color = chain, fill = chain), size = size, alpha = alpha) +
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

#' Criteria to compare Joint Models
#' 
#' Compares two or more fitted joint models using the criteria WAIC, DIC, and LPML.
#' 
#' @param ... two or more objects inheriting from class \code{"jm"}.
#' @param type  which log-likelihood function use to calculate the criteria 
#' (marginal or conditional).  
#' @param order which criteria use to sort the models in the output.
#' @return An object of class \code{compare_jm}. This is a list with the elements: 
#' \tabular{ll}{
#' \code{table} \tab a table with the criteria calculted to each joint model. \cr
#' \code{type} \tab the log-likelihood function used to calculate the criteria.   
#' }
#' @author Dimitris Rizopoulos, \email{d.rizopoulos@@erasmusmc.nl}.
#' @seealso \code{\link{jm}}.
#' @examples 
#' \dontrun{
#' # linear mixed model fits
#' fit_lme1 <- lme(sqrt(CD4) ~ obstime, random = ~ 1 + obstime | patient, 
#'                 data = aids)
#' 
#' fit_lme2 <- lme(sqrt(CD4) ~ obstime + drug, random = ~ 1 + obstime | patient, 
#'                 data = aids)
#' 
#' # cox model fit
#' fit_cox <- coxph(Surv(Time, death) ~ drug, data = aids.id)
#' 
#' # joint model fit 1
#' fit_jm <- jm(fit_cox, fit_lme1, time_var = "obstime")
#' 
#' # joint model fit 2
#' fit_jm2 <- jm(fit_cox, fit_lme2, time_var = "obstime")
#' 
#' # compare the two fitted joint models
#' compare_jm(fit_jm1, fit_jm2)
#' }
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




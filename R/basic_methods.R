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
        for (i in seq_along(parms)) {
            parms_i <- parms[[i]]
            x <- object$mcmc[[parms_i]]
            if (!is.null(x)) coda::gelman.diag(x, ...)
        }
    } else {
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


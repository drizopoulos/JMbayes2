coxph <- function (formula, data, weights, subset, na.action, init, control,
          ties = c("efron", "breslow", "exact"),
          singular.ok = TRUE, robust, model = FALSE, x = FALSE, y = TRUE,
          tt, method = ties, id, cluster, istate, statedata, ...) {
    ties <- match.arg(ties)
    Call <- match.call()
    extraArgs <- list(...)
    if (length(extraArgs)) {
        controlargs <- names(formals(coxph.control))
        indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
        if (any(indx == 0L))
            stop(gettextf("Argument %s not matched", names(extraArgs)[indx ==
                                                                          0L]), domain = NA)
    }
    if (missing(control))
        control <- coxph.control(...)
    if (missing(formula))
        stop("a formula argument is required")
    ss <- c("cluster", "offset")
    if (is.list(formula))
        Terms <- if (missing(data))
            terms(formula[[1]], specials = ss)
    else terms(formula[[1]], specials = ss, data = data)
    else Terms <- if (missing(data))
        terms(formula, specials = ss)
    else terms(formula, specials = ss, data = data)
    tcl <- attr(Terms, "specials")$cluster
    if (length(tcl) > 1)
        stop("a formula cannot have multiple cluster terms")
    if (length(tcl) > 0) {
        factors <- attr(Terms, "factors")
        if (any(factors[tcl, ] > 1))
            stop("cluster() cannot be in an interaction")
        if (attr(Terms, "response") == 0)
            stop("formula must have a Surv response")
        temp <- attr(Terms, "term.labels")
        oo <- attr(Terms, "specials")$offset
        if (!is.null(oo)) {
            ooterm <- rownames(factors)[oo]
            if (oo < tcl)
                temp <- c(ooterm, temp)
            else temp <- c(temp, ooterm)
        }
        if (is.null(Call$cluster))
            Call$cluster <- attr(Terms, "variables")[[1 +
                                                          tcl]][[2]]
        else warning("cluster appears both in a formula and as an argument, formula term ignored")
        #if (is.list(formula))
        #    formula[[1]][[3]] <- reformulate(temp[1 - tcl])[[2]]
        #else formula[[3]] <- reformulate(temp[1 - tcl])[[2]]
        temp. <- temp[1 - tcl]
        if (!length(temp.))
            temp. <- "1"
        if (is.list(formula)) {
            formula[[1]][[3]] <- reformulate(temp.)[[2]]
        } else {
            formula[[3]] <- reformulate(temp.)[[2]]
        }
        Call$formula <- formula
    }
    indx <- match(c("formula", "data", "weights",
                    "subset", "na.action", "cluster", "id",
                    "istate"), names(Call), nomatch = 0)
    if (indx[1] == 0)
        stop("A formula argument is required")
    tform <- Call[c(1, indx)]
    tform[[1L]] <- quote(stats::model.frame)
    if (is.list(formula)) {
        multiform <- TRUE
        dformula <- formula[[1]]
        if (missing(statedata))
            covlist <- parsecovar1(formula[-1])
        else {
            if (!inherits(statedata, "data.frame"))
                stop("statedata must be a data frame")
            if (is.null(statedata$state))
                stop("statedata data frame must contain a 'state' variable")
            covlist <- parsecovar1(formula[-1], names(statedata))
        }
        tlab <- unlist(lapply(covlist$rhs, function(x) attr(terms.formula(x$formula),
                                                            "term.labels")))
        tlab <- c(attr(terms.formula(dformula), "term.labels"),
                  tlab)
        newform <- reformulate(tlab, dformula[[2]])
        environment(newform) <- environment(dformula)
        formula <- newform
        tform$na.action <- na.pass
    }
    else {
        multiform <- FALSE
        covlist <- NULL
        dformula <- formula
    }
    special <- c("strata", "tt")
    tform$formula <- if (missing(data))
        terms(formula, special)
    else terms(formula, special, data = data)
    if (!is.null(attr(tform$formula, "specials")$tt)) {
        coxenv <- new.env(parent = environment(formula))
        assign("tt", function(x) x, envir = coxenv)
        environment(tform$formula) <- coxenv
    }
    mf <- eval(tform, parent.frame())
    if (nrow(mf) == 0)
        stop("No (non-missing) observations")
    Terms <- terms(mf)
    Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv"))
        stop("Response must be a survival object")
    type <- attr(Y, "type")
    multi <- FALSE
    if (type == "mright" || type == "mcounting")
        multi <- TRUE
    else if (type != "right" && type != "counting")
        stop(paste("Cox model doesn't support \"", type,
                   "\" survival data", sep = ""))
    data.n <- nrow(Y)
    if (!multi && multiform)
        stop("formula is a list but the response is not multi-state")
    if (control$timefix)
        Y <- aeqSurv(Y)
    if (length(attr(Terms, "variables")) > 2) {
        ytemp <- terms.inner(formula[1:2])
        xtemp <- terms.inner(formula[-2])
        if (any(!is.na(match(xtemp, ytemp))))
            warning("a variable appears on both the left and right sides of the formula")
    }
    strats <- attr(Terms, "specials")$strata
    hasinteractions <- FALSE
    dropterms <- NULL
    if (length(strats)) {
        stemp <- untangle.specials(Terms, "strata", 1)
        if (length(stemp$vars) == 1)
            strata.keep <- mf[[stemp$vars]]
        else strata.keep <- strata(mf[, stemp$vars], shortlabel = TRUE)
        istrat <- as.integer(strata.keep)
        for (i in stemp$vars) {
            if (any(attr(Terms, "order")[attr(Terms, "factors")[i,
                                                                ] > 0] > 1))
                hasinteractions <- TRUE
        }
        if (!hasinteractions)
            dropterms <- stemp$terms
    }
    else istrat <- NULL
    if (hasinteractions && multi)
        stop("multi-state coxph does not support strata*covariate interactions")
    timetrans <- attr(Terms, "specials")$tt
    if (missing(tt))
        tt <- NULL
    if (length(timetrans)) {
        if (multi)
            stop("the tt() transform is not implemented for multi-state models")
        timetrans <- untangle.specials(Terms, "tt")
        ntrans <- length(timetrans$terms)
        if (is.null(tt)) {
            tt <- function(x, time, riskset, weights) {
                obrien <- function(x) {
                    r <- rank(x)
                    (r - 0.5)/(0.5 + length(r) - r)
                }
                unlist(tapply(x, riskset, obrien))
            }
        }
        if (is.function(tt))
            tt <- list(tt)
        if (is.list(tt)) {
            if (any(!sapply(tt, is.function)))
                stop("The tt argument must contain function or list of functions")
            if (length(tt) != ntrans) {
                if (length(tt) == 1) {
                    temp <- vector("list", ntrans)
                    for (i in 1:ntrans) temp[[i]] <- tt[[1]]
                    tt <- temp
                }
                else stop("Wrong length for tt argument")
            }
        }
        else stop("The tt argument must contain a function or list of functions")
        if (ncol(Y) == 2) {
            if (length(strats) == 0) {
                sorted <- order(-Y[, 1], Y[, 2])
                newstrat <- rep.int(0L, nrow(Y))
                newstrat[1] <- 1L
            }
            else {
                sorted <- order(istrat, -Y[, 1], Y[, 2])
                newstrat <- as.integer(c(1, 1 * (diff(istrat[sorted]) !=
                                                     0)))
            }
            if (storage.mode(Y) != "double")
                storage.mode(Y) <- "double"
            counts <- .Call(Ccoxcount1, Y[sorted, ], as.integer(newstrat))
            tindex <- sorted[counts$index]
        }
        else {
            if (length(strats) == 0) {
                sort.end <- order(-Y[, 2], Y[, 3])
                sort.start <- order(-Y[, 1])
                newstrat <- c(1L, rep(0, nrow(Y) - 1))
            }
            else {
                sort.end <- order(istrat, -Y[, 2], Y[, 3])
                sort.start <- order(istrat, -Y[, 1])
                newstrat <- c(1L, as.integer(diff(istrat[sort.end]) !=
                                                 0))
            }
            if (storage.mode(Y) != "double")
                storage.mode(Y) <- "double"
            counts <- .Call(Ccoxcount2, Y, as.integer(sort.start -
                                                          1L), as.integer(sort.end - 1L), as.integer(newstrat))
            tindex <- counts$index
        }
        Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
        type <- "right"
        mf <- mf[tindex, ]
        istrat <- rep(1:length(counts$nrisk), counts$nrisk)
        weights <- model.weights(mf)
        if (!is.null(weights) && any(!is.finite(weights)))
            stop("weights must be finite")
        tcall <- attr(Terms, "variables")[timetrans$terms +
                                              2]
        pvars <- attr(Terms, "predvars")
        pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
        for (i in 1:ntrans) {
            newtt <- (tt[[i]])(mf[[timetrans$var[i]]], Y[, 1],
                               istrat, weights)
            mf[[timetrans$var[i]]] <- newtt
            nclass <- class(newtt)
            if (any(nclass %in% pmethod)) {
                dummy <- as.call(list(as.name(class(newtt)[1]),
                                      tcall[[i]][[2]]))
                ptemp <- makepredictcall(newtt, dummy)
                pvars[[timetrans$terms[i] + 2]] <- ptemp
            }
        }
        attr(Terms, "predvars") <- pvars
    }
    xlevels <- .getXlevels(Terms, mf)
    cluster <- model.extract(mf, "cluster")
    id <- model.extract(mf, "id")
    if (!is.null(id) && is.null(cluster) && (missing(robust) ||
                                             robust))
        cluster <- id
    if (missing(robust))
        robust <- !is.null(cluster)
    else if (robust && is.null(cluster)) {
        if (ncol(Y) == 2)
            cluster <- seq.int(1, nrow(mf))
        else stop("one of cluster or id is needed")
    }
    contrast.arg <- NULL
    attr(Terms, "intercept") <- 1
    id <- model.extract(mf, "id")
    if (multi) {
        if (length(dropterms)) {
            Terms2 <- Terms[-dropterms]
            dformula <- formula(Terms2)
        }
        else Terms2 <- Terms
        if (length(id) == 0)
            stop("an id statement is required for multi-state models")
        istate <- model.extract(mf, "istate")
        mcheck <- survcheck2(Y, id, istate)
        if (mcheck$flag["overlap"] > 0)
            stop("data set has overlapping intervals for one or more subjects")
        transitions <- mcheck$transitions
        istate <- mcheck$istate
        states <- mcheck$states
        if (missing(statedata))
            covlist2 <- parsecovar2(covlist, NULL, dformula = dformula,
                                    Terms2, transitions, states)
        else covlist2 <- parsecovar2(covlist, statedata, dformula = dformula,
                                     Terms2, transitions, states)
        tmap <- covlist2$tmap
        if (!is.null(covlist)) {
            good.tran <- bad.tran <- rep(FALSE, nrow(Y))
            termname <- rownames(attr(Terms, "factors"))
            trow <- (!is.na(match(rownames(tmap), termname)))
            termiss <- matrix(0L, nrow(mf), ncol(mf))
            for (i in 1:ncol(mf)) {
                xx <- is.na(mf[[i]])
                if (is.matrix(xx))
                    termiss[, i] <- apply(xx, 1, any)
                else termiss[, i] <- xx
            }
            for (i in levels(istate)) {
                rindex <- which(istate == i)
                j <- which(covlist2$mapid[, 1] == match(i, states))
                for (jcol in j) {
                    k <- which(trow & tmap[, jcol] > 0)
                    bad.tran[rindex] <- (bad.tran[rindex] | apply(termiss[rindex,
                                                                          k, drop = FALSE], 1, any))
                    good.tran[rindex] <- (good.tran[rindex] | apply(!termiss[rindex,
                                                                             k, drop = FALSE], 1, all))
                }
            }
            n.partially.used <- sum(good.tran & bad.tran & !is.na(Y))
            omit <- (!good.tran & bad.tran) | is.na(Y)
            if (all(omit))
                stop("all observations deleted due to missing values")
            temp <- setNames(seq(omit)[omit], attr(mf, "row.names")[omit])
            attr(temp, "class") <- "omit"
            mf <- mf[!omit, , drop = FALSE]
            attr(mf, "na.action") <- temp
            Y <- Y[!omit]
            id <- id[!omit]
            if (length(istate))
                istate <- istate[!omit]
        }
    }
    if (length(dropterms)) {
        Terms2 <- Terms[-dropterms]
        X <- model.matrix(Terms2, mf, constrasts = contrast.arg)
        temp <- attr(X, "assign")
        shift <- sort(dropterms)
        temp <- temp + 1 * (shift[1] <= temp)
        if (length(shift) == 2)
            temp + 1 * (shift[2] <= temp)
        attr(X, "assign") <- temp
    }
    else X <- model.matrix(Terms, mf, contrasts = contrast.arg)
    Xatt <- attributes(X)
    if (hasinteractions)
        adrop <- c(0, untangle.specials(Terms, "strata")$terms)
    else adrop <- 0
    xdrop <- Xatt$assign %in% adrop
    X <- X[, !xdrop, drop = FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
    attr(X, "contrasts") <- Xatt$contrasts
    offset <- model.offset(mf)
    if (is.null(offset) | all(offset == 0))
        offset <- rep(0, nrow(mf))
    else if (any(!is.finite(offset) | !is.finite(exp(offset))))
        stop("offsets must lead to a finite risk score")
    weights <- model.weights(mf)
    if (!is.null(weights) && any(!is.finite(weights)))
        stop("weights must be finite")
    assign <- attrassign(X, Terms)
    contr.save <- attr(X, "contrasts")
    if (sum(Y[, ncol(Y)]) == 0) {
        ncoef <- ncol(X)
        ctemp <- rep(NA, ncoef)
        names(ctemp) <- colnames(X)
        concordance = c(concordant = 0, discordant = 0, tied.x = 0,
                        tied.y = 0, tied.xy = 0, concordance = NA, std = NA,
                        timefix = FALSE)
        rval <- list(coefficients = ctemp, var = matrix(0, ncoef,
                                                        ncoef), loglik = c(0, 0), score = 0, iter = 0, linear.predictors = offset,
                     residuals = rep(0, data.n), means = colMeans(X),
                     method = method, n = data.n, nevent = 0, terms = Terms,
                     assign = assign, concordance = concordance, wald.test = 0,
                     y = Y, call = Call)
        class(rval) <- "coxph"
        return(rval)
    }
    if (multi) {
        cmap <- parsecovar3(tmap, colnames(X), attr(X, "assign"))
        xstack <- stacker(cmap, as.integer(istate), X, Y, strata = istrat,
                          states = states)
        X <- xstack$X
        Y <- xstack$Y
        istrat <- xstack$strata
        if (length(offset))
            offset <- offset[xstack$rindex]
        if (length(weights))
            weights <- weights[xstack$rindex]
        if (length(cluster))
            cluster <- cluster[xstack$rindex]
        t2 <- tmap[-1, , drop = FALSE]
        r2 <- row(t2)[!duplicated(as.vector(t2))]
        c2 <- col(t2)[!duplicated(as.vector(t2))]
        a2 <- lapply(seq(along = r2), function(i) {
            cmap[1 + assign[[r2[i]]], c2[i]]
        })
        tab <- table(r2)
        count <- tab[r2]
        names(a2) <- ifelse(count == 1, row.names(t2)[r2], paste(row.names(t2)[r2],
                                                                 colnames(cmap)[c2], sep = "_"))
        assign <- a2
    }
    if (!all(is.finite(X)))
        stop("data contains an infinite predictor")
    if (missing(init))
        init <- NULL
    else {
        if (length(init) != ncol(X))
            stop("wrong length for init argument")
        temp <- X %*% init - sum(colMeans(X) * init) + offset
        if (any(exp(temp) > .Machine$double.xmax) || all(exp(temp) ==
                                                         0))
            stop("initial values lead to overflow or underflow of the exp function")
    }
    pterms <- sapply(mf, inherits, "coxph.penalty")
    if (any(pterms)) {
        pattr <- lapply(mf[pterms], attributes)
        pname <- names(pterms)[pterms]
        ord <- attr(Terms, "order")[match(pname, attr(Terms,
                                                      "term.labels"))]
        if (any(ord > 1))
            stop("Penalty terms cannot be in an interaction")
        pcols <- assign[match(pname, names(assign))]
        fit <- coxpenal.fit(X, Y, istrat, offset, init = init,
                            control, weights = weights, method = method, row.names(mf),
                            pcols, pattr, assign)
    }
    else {
        if (method == "breslow" || method == "efron") {
            if (grepl("right", type))
                fitter <- get("coxph.fit")
            else fitter <- get("agreg.fit")
        }
        else if (method == "exact") {
            if (type == "right")
                fitter <- get("coxexact.fit")
            else fitter <- get("agexact.fit")
        }
        else stop(paste("Unknown method", method))
        fit <- fitter(X, Y, istrat, offset, init, control, weights = weights,
                      method = method, row.names(mf))
    }
    if (is.character(fit)) {
        fit <- list(fail = fit)
        class(fit) <- "coxph"
    }
    else {
        if (!is.null(fit$coefficients) && any(is.na(fit$coefficients))) {
            vars <- (1:length(fit$coefficients))[is.na(fit$coefficients)]
            msg <- paste("X matrix deemed to be singular; variable",
                         paste(vars, collapse = " "))
            if (!singular.ok)
                stop(msg)
        }
        fit$n <- data.n
        fit$nevent <- sum(Y[, ncol(Y)])
        fit$terms <- Terms
        fit$assign <- assign
        class(fit) <- fit$class
        fit$class <- NULL
        if (robust && !is.null(fit$coefficients) && !all(is.na(fit$coefficients))) {
            fit$naive.var <- fit$var
            fit2 <- c(fit, list(x = X, y = Y, weights = weights))
            if (length(istrat))
                fit2$strata <- istrat
            if (length(cluster)) {
                temp <- residuals.coxph(fit2, type = "dfbeta",
                                        collapse = cluster, weighted = TRUE)
                if (is.null(init))
                    fit2$linear.predictors <- 0 * fit$linear.predictors
                else fit2$linear.predictors <- c(X %*% init)
                temp0 <- residuals.coxph(fit2, type = "score",
                                         collapse = cluster, weighted = TRUE)
            }
            else {
                temp <- residuals.coxph(fit2, type = "dfbeta",
                                        weighted = TRUE)
                fit2$linear.predictors <- 0 * fit$linear.predictors
                temp0 <- residuals.coxph(fit2, type = "score",
                                         weighted = TRUE)
            }
            fit$var <- t(temp) %*% temp
            u <- apply(as.matrix(temp0), 2, sum)
            fit$rscore <- coxph.wtest(t(temp0) %*% temp0, u,
                                      control$toler.chol)$test
        }
        if (length(fit$coefficients) && is.null(fit$wald.test)) {
            nabeta <- !is.na(fit$coefficients)
            if (is.null(init))
                temp <- fit$coefficients[nabeta]
            else temp <- (fit$coefficients - init[1:length(fit$coefficients)])[nabeta]
            fit$wald.test <- coxph.wtest(fit$var[nabeta, nabeta],
                                         temp, control$toler.chol)$test
        }
        if (length(cluster))
            temp <- concordancefit(Y, fit$linear.predictors,
                                   istrat, weights, cluster = cluster, reverse = TRUE,
                                   timefix = FALSE)
        else temp <- concordancefit(Y, fit$linear.predictors,
                                    istrat, weights, reverse = TRUE, timefix = FALSE)
        if (is.matrix(temp$count))
            fit$concordance <- c(colSums(temp$count), concordance = temp$concordance,
                                 std = sqrt(temp$var))
        else fit$concordance <- c(temp$count, concordance = temp$concordance,
                                  std = sqrt(temp$var))
        na.action <- attr(mf, "na.action")
        if (length(na.action))
            fit$na.action <- na.action
        if (model) {
            if (length(timetrans)) {
                stop("'model=TRUE' not supported for models with tt terms")
            }
            fit$model <- mf
        }
        if (x) {
            fit$x <- X
            if (length(timetrans))
                fit$strata <- istrat
            else if (length(strats))
                fit$strata <- strata.keep
        }
        if (y)
            fit$y <- Y
        fit$timefix <- control$timefix
    }
    if (!is.null(weights) && any(weights != 1))
        fit$weights <- weights
    names(fit$means) <- names(fit$coefficients)
    if (multi) {
        fit$transitions <- transitions
        fit$states <- states
        fit$cmap <- cmap
        fit$resid <- rowsum(fit$resid, xstack$rindex)
        names(fit$coefficients) <- seq(along = fit$coefficients)
        if (x)
            fit$strata <- istrat
        class(fit) <- c("coxphms", class(fit))
    }
    fit$formula <- formula(Terms)
    if (length(xlevels) > 0)
        fit$xlevels <- xlevels
    fit$contrasts <- contr.save
    if (any(offset != 0))
        fit$offset <- offset
    fit$call <- Call
    fit
}

extract_terms <- function (object, which = c("fixed", "random"), data) {
    which <- match.arg(which)
    if (inherits(object, "MixMod")) {
        if (which == "fixed") object$Terms$termsX else object$Terms$termsZ[[1L]]
    } else {
        if (which == "fixed") object$terms else {
            # here we extract first grouping variable;
            # it won't work for nested random effects
            form <- attr(object$modelStruct$reStruct[[1L]], "formula")
            terms(model.frame(form, data = data))
        }
    }
}

fix_NAs_fixed <- function (x, NAs_fixed, NAs_random) {
    if (is.null(NAs_random)) {
        x
    } else {
        extra_nas <- NAs_random[!NAs_random %in% NAs_fixed]
        if (length(extra_nas)) {
            if (is.data.frame(x)) x[-extra_nas, , drop = FALSE] else x[-extra_nas]
        } else x
    }
}

fix_NAs_random <- function (z, NAs_random, NAs_fixed) {
    if (is.null(NAs_fixed)) {
        z
    } else {
        extra_nas <- NAs_fixed[!NAs_fixed %in% NAs_random]
        if (length(extra_nas)) {
            if (is.data.frame(z)) z[-extra_nas, , drop = FALSE] else z[-extra_nas]
        } else z
    }
}

exclude_NAs <- function (NAs_FE, NAs_RE, id) {
    all_NAs <- union(NAs_FE, NAs_RE)
    if (!is.null(all_NAs)) id[-all_NAs] else id
}

bdiag <- function (...) {
    mlist <- list(...)
    if (length(mlist) == 1)
        mlist <- unlist(mlist, recursive = FALSE)
    csdim <- rbind(c(0, 0), apply(sapply(mlist, dim), 1, cumsum))
    ret <- array(0, dim = csdim[length(mlist) + 1, ])
    add1 <- matrix(rep(1:0, 2), ncol = 2)
    for (i in seq_along(mlist)) {
        indx <- apply(csdim[i:(i + 1), ] + add1, 2, function(x) x[1]:x[2])
        if (is.null(dim(indx))) {
            ret[indx[[1]], indx[[2]]] <- mlist[[i]]
        }
        else {
            ret[indx[, 1], indx[, 2]] <- mlist[[i]]
        }
    }
    colnames(ret) <- unlist(lapply(mlist, colnames))
    ret
}

right_rows <- function (data, times, ids, Q_points) {
    fids <- factor(ids, levels = unique(ids))
    if (!is.list(Q_points))
        Q_points <- split(Q_points, row(Q_points))
    ind <- mapply(findInterval, Q_points, split(times, fids))
    ind[ind < 1] <- 1
    rownams_id <- split(row.names(data), fids)
    ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
    data[c(ind), ]
}

extract_functional_forms_per_outcome <- function (i) {
    Form <- formula(functional_form, rhs = i)
    term_labels <- attr(terms(Form), "term.labels")
    possible_forms <- c("value", "slope", "area")
    sapply(possible_forms, grep, x = term_labels, fixed = TRUE)
}

LongData_HazardModel <- function (time_points, data, times, ids, timeVar) {
    unq_ids <- unique(ids)
    fids <- factor(ids, levels = unq_ids)
    if (!is.matrix(time_points)) {
        time_points <- as.matrix(time_points)
    }
    if (nrow(time_points) != length(unq_ids)) {
        stop("the length of unique 'ids' does not match the number of rows ",
             "of 'time_points'.")
    }
    tt <- split(time_points, row(time_points))
    ind <- mapply(findInterval, tt, split(times, fids))
    ind[ind < 1] <- 1
    if (!is.matrix(ind)) {
        ind <- rbind(ind)
    }
    rownams_id <- split(row.names(data), fids)
    ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
    data <- data[c(ind), ]
    data[[timeVar]] <- c(t(time_points))
    row.names(data) <- seq_len(nrow(data))
    data
}

last_rows <- function (data, ids) {
    fidVar <- factor(ids, levels = unique(ids))
    data[tapply(row.names(data), fidVar, tail, n = 1L), ]
}

gaussKronrod <- function (k = 15L) {
    sk <- c(-0.949107912342758, -0.741531185599394, -0.405845151377397,
            0, 0.405845151377397, 0.741531185599394, 0.949107912342758,
            -0.991455371120813, -0.864864423359769, -0.586087235467691,
            -0.207784955007898, 0.207784955007898, 0.586087235467691,
            0.864864423359769, 0.991455371120813)
    wk15 <- c(0.0630920926299786, 0.140653259715526, 0.190350578064785,
              0.209482141084728, 0.190350578064785, 0.140653259715526,
              0.0630920926299786, 0.0229353220105292, 0.10479001032225,
              0.169004726639268, 0.204432940075299, 0.204432940075299,
              0.169004726639268, 0.10479001032225, 0.0229353220105292)
    wk7 <- c(0.12948496616887, 0.279705391489277, 0.381830050505119,
             0.417959183673469, 0.381830050505119, 0.279705391489277,
             0.12948496616887)
    if (k == 7L) {
        list(sk = sk[seq_len(k)], wk = wk7)
    } else {
        list(sk = sk, wk = wk15)
    }
}

desing_matrices_functional_forms <- function (time, terms, data, timeVar, idVar,
                                              Fun_Forms) {
    desgn_matr <- function (time, terms) {
        D <- LongData_HazardModel(time, data, data[[timeVar]],
                                  data[[idVar]], timeVar)
        mf <- lapply(terms, model.frame.default, data = D)
        mapply(model.matrix.default, terms, mf)
    }
    degn_matr_slp <- function (time, terms) {
        M1 <- desgn_matr(time + 0.001, terms)
        M2 <- desgn_matr(time - 0.001, terms)
        mapply(function (x1, x2) (x1 - x2) / 0.002, M1, M2)
    }
    degn_matr_area <- function (time, terms) {
        if (is.matrix(time)) {
            time <- c(t(time))
        }
        GK <- gaussKronrod(15L)
        wk <- GK$wk
        sk <- GK$sk
        P <- unname(time / 2)
        st <- outer(P, sk + 1)
        out <- vector("list", 15L)
        for (i in seq_len(15L)) {
            ss <- if (nrow(st) == length(unique(data[[idVar]]))) {
                st[, i]
            } else {
                matrix(st[, i], ncol = 15, byrow = TRUE)
            }
            M <- desgn_matr(ss, terms)
            out[[i]] <- lapply(M, "*", P * wk[i])
        }
        lapply(seq_along(M), function (i) Reduce("+", lapply(out, "[[", i)))
    }
    ################
    out <- list("value" = desgn_matr(time, terms),
                "slope" = degn_matr_slp(time, terms),
                "area" = degn_matr_area(time, terms))
    out <- lapply(seq_along(Fun_Forms), function (i) lapply(out[Fun_Forms[[i]]], "[[", i))
    names(out) <- names(Fun_Forms)
    out
}

##########################################################################################
# GP: ADDED ALTERNATIVE VERSION OF desing_matrices_functional_forms() THAT RETURNS ARRAY #
# INSTEAD OF LIST ########################################################################

desing_matrices_functional_forms_array <- function (time, terms, data, timeVar, idVar,
                                                    Fun_Forms) {
    desgn_matr <- function (time, terms) {
        D <- LongData_HazardModel(time, data, data[[timeVar]],
                                  data[[idVar]], timeVar)
        mf <- lapply(terms, model.frame.default, data = D)
        mapply(model.matrix.default, terms, mf)
        }
    degn_matr_slp <- function (time, terms) {
        M1 <- desgn_matr(time + 0.001, terms)
        M2 <- desgn_matr(time - 0.001, terms)
        mapply(function (x1, x2) (x1 - x2) / 0.002, M1, M2)
        }
    degn_matr_area <- function (time, terms) {
        if (is.matrix(time)) {
            time <- c(t(time))
        }
        GK <- gaussKronrod(15L)
        wk <- GK$wk
        sk <- GK$sk
        P <- unname(time / 2)
        st <- outer(P, sk + 1)
        out <- vector("list", 15L)
        for (i in seq_len(15L)) {
            ss <- if (nrow(st) == length(unique(data[[idVar]]))) {
                st[, i]
            } else {
                matrix(st[, i], ncol = 15, byrow = TRUE)
            }
            M <- desgn_matr(ss, terms)
            out[[i]] <- lapply(M, "*", P * wk[i])
        }
        lapply(seq_along(M), function (i) Reduce("+", lapply(out, "[[", i)))
    }
    ################
    out <- list("value" = desgn_matr(time, terms),
                "slope" = degn_matr_slp(time, terms),
                "area" = degn_matr_area(time, terms))
    out <- lapply(seq_along(Fun_Forms), function (i) lapply(out[Fun_Forms[[i]]], "[[", i))
    names(out) <- names(Fun_Forms)
    lapply(out, function (x) array(unlist(x), dim = c(nrow(x[[1]]), ncol(x[[1]]), length(x))))
}

extract_D <- function (object) {
    if (inherits(object, "lme")) {
        lapply(pdMatrix(object$modelStruct$reStruct), "*",
               object$sigma^2)[[1]]
    } else {
        object$D
    }
}

knots <- function (xl, xr, ndx, deg) {
    # From Paul Eilers
    dx <- (xr - xl) / ndx
    seq(xl - deg * dx, xr + deg * dx, by = dx)
}

SurvData_HazardModel <- function (time_points, data, times, ids) {
    unq_ids <- unique(ids)
    fids <- factor(ids, levels = unq_ids)
    if (!is.matrix(time_points)) {
        time_points <- as.matrix(time_points)
    }
    if (nrow(time_points) != length(unq_ids)) {
        stop("the length of unique 'ids' does not match the number of rows ",
             "of 'time_points'.")
    }
    tt <- split(time_points, row(time_points))
    ind <- mapply(findInterval, tt, split(times, fids))
    ind[ind < 1] <- 1
    if (!is.matrix(ind)) {
        ind <- rbind(ind)
    }
    rownams_id <- split(row.names(data), fids)
    ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
    data <- data[c(ind), ]
}

extract_b <- function (object, id, n) {
    b <- data.matrix(ranef(object))
    mat <- matrix(0.0, n, ncol(b))
    colnames(mat) <- colnames(b)
    mat[id, ] <- b
    mat
}

extract_log_sigmas <- function (object) {
    if (inherits(object, "lme")) {
        # we extract the log of sigma to be consisten with GLMMadaptive
        log(object$sigma)
    } else {
        object$phis
    }
}

value <- slope <- area <- function (x) rep(1, length(x))


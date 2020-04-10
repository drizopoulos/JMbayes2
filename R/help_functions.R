cd <- function (x, f, ..., eps = 0.001) {
    # numerical derivative: central difference approximation for scalar functions
    n <- length(x)
    res <- numeric(n)
    ex <- pmax(abs(x), 1)
    for (i in seq_len(n)) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[i] <- diff.f/diff.x
    }
    res
}

cd_vec <- function (x, f, ..., eps = 0.001) {
    # numerical derivative: central difference approximation for vector functions
    n <- length(x)
    res <- matrix(0, n, n)
    ex <- pmax(abs(x), 1)
    for (i in seq_len(n)) {
        x1 <- x2 <- x
        x1[i] <- x[i] + eps * ex[i]
        x2[i] <- x[i] - eps * ex[i]
        diff.f <- c(f(x1, ...) - f(x2, ...))
        diff.x <- x1[i] - x2[i]
        res[, i] <- diff.f/diff.x
    }
    0.5 * (res + t(res))
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
    # constructs a block-diagonal matrix
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

.bdiag <- function (mlist) {
    # constructs a block-diagonal matrix
    mlist <- mlist[sapply(mlist, length) > 0]
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
    # find the right row to use from a data.frame. This is used to find which
    # specific rows from the longitudinal or survival datasets need to be used in
    # the specification of the design matrices for the survival model.
    fids <- factor(ids, levels = unique(ids))
    if (!is.list(Q_points))
        Q_points <- split(Q_points, row(Q_points))
    ind <- mapply(findInterval, Q_points, split(times, fids))
    ind[ind < 1] <- 1
    rownams_id <- split(row.names(data), fids)
    ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
    data[c(ind), ]
}

extract_functional_forms_per_outcome <- function (Form) {
    term_labels <- attr(terms(Form), "term.labels")
    possible_forms <- c("value", "slope", "area")
    sapply(possible_forms, grep, x = term_labels, fixed = TRUE)
}

LongData_HazardModel <- function (time_points, data, times, ids, timeVar) {
    unq_ids <- unique(ids)
    fids <- factor(ids, levels = unq_ids)
    tt <- if (is.list(time_points)) {
        time_points
    } else {
        if (!is.matrix(time_points)) {
            time_points <- as.matrix(time_points)
        }
        split(time_points, row(time_points))
    }
    if (length(tt) != length(unq_ids)) {
        stop("the length of unique 'ids' does not match the number of rows ",
             "of 'time_points'.")
    }
    ind <- mapply(findInterval, tt, split(times, fids))
    rownams_id <- split(row.names(data), fids)
    if (!is.list(ind)) {
        ind[ind < 1] <- 1
        if (!is.matrix(ind)) {
            ind <- rbind(ind)
        }
        ind <- mapply(`[`, rownams_id, split(ind, col(ind)), SIMPLIFY = FALSE)
    } else {
        ind <- lapply(ind, function (x) {x[x < 1] <- 1; x})
        ind <- mapply(`[`, rownams_id, ind, SIMPLIFY = FALSE)
    }
    data <- data[unlist(ind, use.names = FALSE), ]
    data[[timeVar]] <- if (is.matrix(time_points)) {
        c(t(time_points))
    } else {
        unlist(time_points, use.names = FALSE)
    }
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
        if (is.list(time)) {
            t1 <- lapply(time, function (t) t + 0.001)
            t2 <- lapply(time, function (t) t - 0.001)
            M1 <- desgn_matr(t1, terms)
            M2 <- desgn_matr(t2, terms)
        } else {
            M1 <- desgn_matr(time + 0.001, terms)
            M2 <- desgn_matr(time - 0.001, terms)
        }
        mapply(function (x1, x2) (x1 - x2) / 0.002, M1, M2)
    }
    degn_matr_area <- function (time, terms) {
        if (!is.list(time)) {
            time <- if (is.matrix(time)) split(time, row(time))
            else split(time, seq_along(time))
        }
        GK <- gaussKronrod(15L)
        wk <- GK$wk
        sk <- GK$sk
        quadrature_points <- function (x) {
            P <- unname(x / 2)
            sk <- outer(P, sk + 1)
            list(P = c(t(outer(P, wk))), sk = sk)
        }
        qp <- lapply(time, quadrature_points)
        ss <- lapply(qp, function (x) c(t(x[['sk']])))
        Pwk <- unlist(lapply(qp, '[[', 'P'), use.names = FALSE)
        M <- desgn_matr(ss, terms)
        M <- lapply(M, "*", Pwk)
        sum_qp <- function (m) {
            n <- nrow(m)
            grp <- rep(seq_len(round(n / 15)), each = 15L)
            rowsum(m, grp, reorder = FALSE)
        }
        lapply(M, sum_qp)
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
    tt <- if (is.list(time_points)) {
        time_points
    } else {
        if (!is.matrix(time_points)) {
            time_points <- as.matrix(time_points)
        }
        split(time_points, row(time_points))
    }
    if (length(tt) != length(unq_ids)) {
        stop("the length of unique 'ids' does not match the number of rows ",
             "of 'time_points'.")
    }
    ind <- mapply(findInterval, tt, split(times, fids))
    rownams_id <- split(row.names(data), fids)
    if (!is.list(ind)) {
        ind[ind < 1] <- 1
        if (!is.matrix(ind)) {
            ind <- rbind(ind)
        }
        ind <- mapply(`[`, rownams_id, split(ind, col(ind)), SIMPLIFY = FALSE)
    } else {
        ind <- lapply(ind, function (x) {x[x < 1] <- 1; x})
        ind <- mapply(`[`, rownams_id, ind, SIMPLIFY = FALSE)
    }
    data <- data[unlist(ind, use.names = FALSE), ]
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
        # we extract the log of sigma to be consistent with GLMMadaptive
        log(object$sigma)
    } else {
        object$phis
    }
}

value <- slope <- area <- function (x) rep(1, length(x))

create_HC_X <- function (TermsX, TermsZ, x, z, id, mfHC) {
    # function that creates the hierarchical centering version of the
    # design matrix for the fixed effects
    find_positions <- function (nams1, nams2) {
        nams1 <- gsub("^", "\\^", nams1, fixed = TRUE)
        vals <- c(glob2rx(nams1), glob2rx(paste0(nams1, ":*")),
                  glob2rx(paste0("*:", nams1)))
        out <- sort(unique(unlist(lapply(vals, grep, x = nams2))))
        out
    }
    check_td <- function (x, id) {
        !all(sapply(split(x, id), function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
    }
    grep2 <- function (x, nams) grep(x, nams, fixed = TRUE)
    terms.labs_X <- attr(TermsX, "term.labels")
    terms.labs_Z <- attr(TermsZ, "term.labels")
    # check for time-varying covariates
    timeTerms <- if (length(terms.labs_Z)) {
        unlist(lapply(terms.labs_Z, grep2, nams = colnames(x)))
    }
    which_td <- unname(which(apply(x, 2, check_td, id = id)))
    all_TDterms <- unique(c(timeTerms, which_td))
    baseline <- if (length(all_TDterms)) seq_len(ncol(x))[-all_TDterms] else seq_len(ncol(x))
    ind_colmns <- c(list(baseline), lapply(colnames(z)[-1L], find_positions,
                                           nams2 = colnames(x)))
    ind_colmns2 <- seq_len(ncol(x))
    ind_colmns2 <- ind_colmns2[!ind_colmns2 %in% unlist(ind_colmns)]
    mfHC <- mfHC[!duplicated(id), ]
    Xhc <- if (length(terms.labs_Z)) {
        which.timevar <- unique(unlist(lapply(terms.labs_Z, grep2, nams = names(mfHC))))
        mfHC[which.timevar] <- lapply(mfHC[which.timevar],
                                      function (x) { x[] <- 1; x })
        model.matrix(TermsX, mfHC)
    } else {
        model.matrix(TermsX, mfHC)
    }
    index <- numeric(ncol(Xhc))
    for (i in seq_along(ind_colmns)) {
        index[ind_colmns[[i]]] <- i
    }
    list(Xhc = Xhc, columns_HC = index, columns_nHC = ind_colmns2)
}

drop_names <- function (x) {
    unname2 <- function (x) {
        if (!is.null(attr(x, "assign"))) {
            attr(x, "assign") <- NULL
        }
        if (!is.null(attr(x, "contrasts"))) {
            attr(x, "contrasts") <- NULL
        }
        unname(x)
    }
    drp <- function (z) if (is.list(z)) lapply(z, unname2) else unname2(z)
    if (is.list(x)) lapply(x, drp) else unname2(x)
}

vcov2 <- function (model) {
    if (inherits(model, "MixMod")) vcov(model, "fixed-effects") else vcov(model)
}

get_vcov_FE <- function (model, cc, which = c("betas", "tilde_betas")) {
    ind <- cc
    if (which == "betas") {
        return(if (length(ind)) vcov2(model)[-ind, -ind, drop = FALSE] else vcov2(model))
    }
    if (which == "tilde_betas" && length(ind)) vcov2(model)[ind, ind, drop = FALSE] else NULL
}

extract_vcov_prop_RE <- function (object, Z_k, id_k) {
    # variance-covariance matrix for the proposal distributions for the random effects
    # (one proposal distribution per subject/cluster)
    if (inherits(object, "lme")) {
        D <- lapply(pdMatrix(object$modelStruct$reStruct), "*",
                    object$sigma^2)[[1]]
        invD <- solve(D)
        sigma <- object$sigma
        sigma2 <- sigma * sigma
        n <- length(unique(id_k))
        cov_postRE <- vector("list", n)
        names(cov_postRE) <- unique(id_k)
        for (i in seq_len(n)) {
            Z_k_i <- Z_k[id_k == i, , drop = FALSE]
            cov_postRE[[i]] <- solve.default(crossprod(Z_k_i) / sigma2 + invD)
        }
        cov_postRE
    } else if (inherits(object, "MixMod")) {
        out <- object$post_vars
        names(out) <- unique(id_k)
        out
    }
}

init_vals_surv <- function(Data, model_info, data, betas, b, control) {
    # initial values and variance-covariance matrix for the proposal distributions
    # for the coefficients in the survival model. The function does the following two
    # things: (1) Fits a time-varying Cox model using the baseline covariates from
    # 'Surv_object' and the longitudinal outcomes from 'Mixed_objects', taking also into
    # account the functional forms. From this first step we get initial values the
    # variance-covariance matrix for the proposals in the MCMC for the 'gammas' and
    # 'alphas' coefficients. Then we go in step (2): We specify a function that calculates
    # the log likelihood for 'bs_gammas' (i.e., the log density of the survival submodel).
    # In this step, and contrary to the previous step we do account for left or interval
    # censoring. We get the initial values for 'bs_gammas' and the corresponding vcov
    # matrix for the proposal in the MCMC using optim().
    Time_start <- Data$Time_start
    Time_right <- Data$Time_right
    Time_left <- Data$Time_left
    delta <- Data$delta
    which_event <- Data$which_event
    which_right <- Data$which_right
    which_left <- Data$which_left
    which_interval <- Data$which_interval
    if (length(which_interval)) {
        Time_right[which_interval] <- 0.5 * (Time_right[which_interval] +
                                                 Time_left[which_interval])
        delta <- as.numeric(delta == 1 | delta == 3)
    }
    if (length(which_left)) {
        Time_right[which_left] <- Time_left[which_left]
    }
    n <- model_info$n
    ###
    dataL <- data$dataL
    dataS <- data$dataS
    ###
    idVar <- model_info$var_names$idVar
    time_var <- model_info$var_names$time_var
    idT <- model_info$ids$idT
    terms_FE_noResp <- model_info$terms$terms_FE_noResp
    terms_RE <- model_info$terms$terms_RE
    terms_RE <- model_info$terms$terms_RE
    terms_Surv_noResp <- model_info$terms$terms_Surv_noResp
    ###
    functional_forms <- model_info$fun_forms$functional_forms
    functional_forms_per_outcome <- model_info$fun_forms$functional_forms_per_outcome
    collapsed_functional_forms <- model_info$fun_forms$collapsed_functional_forms
    ###
    X_H <- Data$X_H; Z_H <- Data$Z_H; U_H <- Data$U_H; W0_H <- Data$W0_H; W_H <- Data$W_H
    X_h <- Data$X_h; Z_h <- Data$Z_h; U_h <- Data$U_h; W0_h <- Data$W0_h; W_h <- Data$W_h
    X_H2 <- Data$X_H2; Z_H2 <- Data$Z_H2; U_H2 <- Data$U_H2; W0_H2 <- Data$W0_H2; W_H2 <- Data$W_H2
    ###
    log_Pwk <- Data[["log_Pwk"]]
    log_Pwk2 <- Data[["log_Pwk2"]]
    ######################################################################################
    ######################################################################################
    times_long <- split(dataL[[time_var]], dataL[[idVar]])
    dataS_init <- SurvData_HazardModel(times_long, dataS, Time_start, idT)
    mf <- model.frame.default(terms_Surv_noResp, data = dataS_init)
    W_init <- model.matrix.default(terms_Surv_noResp, mf)[, -1, drop = FALSE]
    if (!ncol(W_init)) {
        W_init <- cbind(W_init, rep(0, nrow(W_init)))
    }
    X_init <- desing_matrices_functional_forms(times_long, terms_FE_noResp,
                                               dataL, time_var, idVar,
                                               collapsed_functional_forms)
    Z_init <- desing_matrices_functional_forms(times_long, terms_RE,
                                               dataL, time_var, idVar,
                                               collapsed_functional_forms)
    U_init <- lapply(seq_along(X_init), function (i) {
        tt <- terms(functional_forms[[i]])
        model.matrix(tt, model.frame(tt, data = dataS_init))[, -1, drop = FALSE]
    })
    ##############
    id_init <- rep(list(dataL[[idVar]]), length.out = length(X_init))
    eta_init <- linpred_surv(X_init, betas, Z_init, b, id_init)
    Wlong_init <- create_Wlong(eta_init, functional_forms_per_outcome, U_init)
    Wlong_init <- do.call("cbind", Wlong_init)
    ######################################################################################
    ######################################################################################
    start <- dataL[[time_var]]
    fid <- dataL[[idVar]]
    fid <- factor(fid, levels = unique(fid))
    stop <- unlist(mapply(`c`, tapply(start, fid, tail, n = -1), split(Time_right, idT),
                          SIMPLIFY = FALSE), use.names = FALSE)
    create_event <- function (ni, delta) {
        if (ni == 1) delta else c(rep(0, ni - 1), delta)
    }
    event <- unlist(mapply(create_event, ni = tapply(fid, fid, length), delta,
                           SIMPLIFY = FALSE), use.names = FALSE)
    any_gammas <- !(ncol(W_init) == 1 && all(W_init[, 1] == 0))
    WW <- if (any_gammas) cbind(W_init, Wlong_init) else Wlong_init
    ####
    fm <- coxph(Surv(start, stop, event) ~ WW)
    coefs <- coef(fm)
    gammas <- if (any_gammas) head(coefs, ncol(W_init)) else 0.0
    alphas <- tail(coefs, ncol(Wlong_init))
    alphas <- split(alphas, rep(seq_along(U_H), sapply(U_H, ncol)))
    V <- vcov(fm)
    if (any_gammas) {
        vcov_prop_gammas <- V[1:ncol(W_init), 1:ncol(W_init)]
        vcov_prop_alphas <- V[-(1:ncol(W_init)), -(1:ncol(W_init))]
    } else {
        vcov_prop_gammas <- matrix(0.0, 1, 1)
        vcov_prop_alphas <- V
    }
    out <- list(gammas = gammas, alphas = alphas,
                vcov_prop_gammas = vcov_prop_gammas, vcov_prop_alphas = vcov_prop_alphas)
    ######################################################################################
    ######################################################################################
    id_H <- lapply(X_H, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
    eta_H <- linpred_surv(X_H, betas, Z_H, b, id_H)
    Wlong_H <- create_Wlong(eta_H, functional_forms_per_outcome, U_H)
    if (length(which_event)) {
        id_h <- lapply(X_h, function (x) seq_len(nrow(x[[1]])))
        eta_h <- linpred_surv(X_h, betas, Z_h, b, id_h)
        Wlong_h <- create_Wlong(eta_h, functional_forms_per_outcome, U_h)
    } else {
        Wlong_h <- rep(list(matrix(0.0, length(Time_right), 1)), length(W_H))
    }
    if (length(which_interval)) {
        id_H2 <- lapply(X_H2, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
        eta_H2 <- linpred_surv(X_H2, betas, Z_H, b, id_H2)
        Wlong_H2 <- create_Wlong(eta_H2, functional_forms_per_outcome, U_H2)
    } else {
        Wlong_H2 <- rep(list(matrix(0.0, length(Time_right), 1)), length(W_H))
    }

    log_dens_surv <- function (bs_gammas) {
        lambda_H <- W0_H %*% bs_gammas + W_H %*% gammas
        for (i in seq_along(Wlong_H)) {
            lambda_H <- lambda_H + Wlong_H[[i]] %*% alphas[[i]]
        }
        lambda_h <- matrix(0.0, n, 1)
        if (length(which_event)) {
            lambda_h <- W0_h %*% bs_gammas + W_h %*% gammas
            for (i in seq_along(Wlong_h)) {
                W_h_i <- Wlong_h[[i]]
                lambda_h <- lambda_h + W_h_i %*% alphas[[i]]
            }
        }
        lambda_H2 <- matrix(0.0, nrow(Wlong_H2[[1]]), 1)
        if (length(which_interval)) {
            lambda_H2 <- W0_H2 %*% bs_gammas + W_H2 %*% gammas
            for (i in seq_along(Wlong_H2)) {
                W_H2_i <- Wlong_H2[[i]]
                lambda_H2 <- lambda_H2 + W_H2_i %*% alphas[[i]]
            }
        }
        H <- rowsum(exp(log_Pwk + lambda_H), group = id_H[[1]], reorder = FALSE)
        log_Lik_surv <- numeric(n)
        which_right_event <- c(which_right, which_event)
        if (length(which_right_event)) {
            log_Lik_surv[which_right_event] <- - H[which_right_event]
        }
        if (length(which_event)) {
            log_Lik_surv[which_event] <- log_Lik_surv[which_event] + lambda_h[which_event]
        }
        if (length(which_left)) {
            log_Lik_surv[which_left] <- log1p(- exp(- H[which_left]))
        }
        if (length(which_interval)) {
            H2 <- rowsum(exp(log_Pwk2 + lambda_H2), group = id_H2[[1]], reorder = FALSE)
            log_Lik_surv[which_interval] <- log(exp(- H[which_interval]) -
                                                    exp(- (H2[which_interval] + H[which_interval])))
        }
        - sum(log_Lik_surv, na.rm = TRUE)
    }
    opt <- optim(rep(0.01, ncol(W0_H)), log_dens_surv, method = "BFGS",
                 hessian = TRUE)
    c(out, list(bs_gammas = opt$par, vcov_prop_bs_gammas = solve(opt$hessian)))
}

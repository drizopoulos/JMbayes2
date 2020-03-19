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

.bdiag <- function (mlist) {
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


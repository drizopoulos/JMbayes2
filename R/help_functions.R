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

symm_mat <- function (M) {
    0.5 * (M + t(M))
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
    #if (length(mlist) == 1)
    #    mlist <- unlist(mlist, recursive = FALSE)
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
    #colnames(ret) <- unlist(lapply(mlist, colnames))
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

extract_functional_forms <- function (Form, data) {
    tr <- terms(Form)
    mF <- model.frame(tr, data = data)
    M <- model.matrix(tr, mF)
    cnams <- colnames(M)
    possible_forms <- c("value(", "slope(", "area(")
    ind <- unlist(lapply(possible_forms, grep, x = cnams, fixed = TRUE))
    M <- M[, cnams %in% cnams[unique(ind)], drop = FALSE]
    sapply(c("value", "slope", "area"), grep, x = colnames(M), fixed = TRUE)
}

construct_Umat <- function (fForms, dataS) {
    tt <- terms(fForms)
    m <- model.matrix(tt, model.frame(tt, data = dataS))
    cnams <- colnames(m)
    ind_value <- grep("value(", cnams, fixed = TRUE)
    ind_slope <- grep("slope(", cnams, fixed = TRUE)
    ind_area <- grep("area(", cnams, fixed = TRUE)
    ind <- unique(c(ind_value, ind_slope, ind_area))
    m[, cnams %in% cnams[ind], drop = FALSE]
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
    #if (length(tt) != length(unq_ids)) {
    #    stop("the length of unique 'ids' does not match the number of rows ",
    #         "of 'time_points'.")
    #}
    ind <- mapply(findInterval, tt, split(times, fids))
    rownams_id <- split(row.names(data), fids)
    if (!is.list(ind)) {
        ind[ind < 1] <- 1
        if (!is.matrix(ind)) {
            ind <- rbind(ind)
        }
        ind <- mapply2(`[`, rownams_id, split(ind, col(ind)))
    } else {
        ind <- lapply(ind, function (x) {x[x < 1] <- 1; x})
        ind <- mapply2(`[`, rownams_id, ind)
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

locf <- function (object, fromLast = FALSE, maxgap = Inf) {
    .fill_short_gaps <- function (x, fill, maxgap) {
        if (maxgap <= 0)
            return(x)
        if (maxgap >= length(x))
            return(fill)
        naruns <- rle(is.na(x))
        naruns$values[naruns$lengths > maxgap] <- FALSE
        naok <- inverse.rle(naruns)
        x[naok] <- fill[naok]
        return(x)
    }
    if (fromLast) object <- rev(object)
    ok <- which(!is.na(object))
    if (is.na(object[1L]))
        ok <- c(1L, ok)
    gaps <- diff(c(ok, length(object) + 1L))
    object <- if (any(gaps > maxgap)) {
        .fill_short_gaps(object, rep(object[ok], gaps), maxgap = maxgap)
    } else {
        rep(object[ok], gaps)
    }
    if (fromLast)
        object <- rev(object)
    object
}

design_matrices_functional_forms <- function (time, terms, data, timeVar, idVar,
                                              Fun_Forms, Xbar = NULL) {
    data[] <- lapply(data, function (x) locf(locf(x), fromLast = TRUE))
    desgn_matr <- function (time, terms, Xbar) {
        D <- LongData_HazardModel(time, data, data[[timeVar]],
                                  data[[idVar]], timeVar)
        mf <- lapply(terms, model.frame.default, data = D)
        X <- mapply2(model.matrix.default, terms, mf)
        if (!is.null(Xbar))
            X <- mapply2(function (m, mu) m - rep(mu, each = nrow(m)), X, Xbar)
        X
    }
    degn_matr_slp <- function (time, terms, Xbar) {
        if (is.list(time)) {
            t1 <- lapply(time, function (t) t + 0.001)
            t2 <- lapply(time, function (t) t - 0.001)
            M1 <- desgn_matr(t1, terms, Xbar)
            M2 <- desgn_matr(t2, terms, Xbar)
        } else {
            M1 <- desgn_matr(time + 0.001, terms, Xbar)
            M2 <- desgn_matr(time - 0.001, terms, Xbar)
        }
        mapply2(function (x1, x2) (x1 - x2) / 0.002, M1, M2)
    }
    degn_matr_area <- function (time, terms, Xbar) {
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
            # we divide with x to obtain the standardized area
            list(P = c(t(outer(P / x, wk))), sk = sk)
        }
        qp <- lapply(time, quadrature_points)
        ss <- lapply(qp, function (x) c(t(x[['sk']])))
        Pwk <- unlist(lapply(qp, '[[', 'P'), use.names = FALSE)
        M <- desgn_matr(ss, terms, Xbar)
        M <- lapply(M, "*", Pwk)
        sum_qp <- function (m) {
            n <- nrow(m)
            grp <- rep(seq_len(round(n / 15)), each = 15L)
            rowsum(m, grp, reorder = FALSE)
        }
        lapply(M, sum_qp)
    }
    ################
    out <- list("value" = desgn_matr(time, terms, Xbar),
                "slope" = degn_matr_slp(time, terms, Xbar),
                "area" = degn_matr_area(time, terms, Xbar))
    out <- lapply(seq_along(Fun_Forms), function (i) lapply(out[Fun_Forms[[i]]], "[[", i))
    names(out) <- names(Fun_Forms)
    out
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

SurvData_HazardModel <- function (time_points, data, times, ids, time_var) {
    unq_ids <- unique(ids)
    fids <- factor(ids, levels = unq_ids)
    tt <- if (is.list(time_points)) {
        time_points
    } else {
        if (!is.matrix(time_points)) {
            time_points <- as.matrix(time_points)
        }
        sp <- split(time_points, row(time_points))
        names(sp) <- rownames(time_points)
        sp
    }
    spl_times <- split(times, fids)[names(tt)]
    ind <- mapply(findInterval, tt, spl_times)
    rownams_id <- split(row.names(data), fids)
    if (!is.list(ind)) {
        ind[ind < 1] <- 1
        if (!is.matrix(ind)) {
            ind <- rbind(ind)
        }
        ind <- mapply2(`[`, rownams_id, split(ind, col(ind)))
    } else {
        ind <- lapply(ind, function (x) {x[x < 1] <- 1; x})
        ind <- mapply2(`[`, rownams_id, ind)
    }
    data <- data[unlist(ind, use.names = FALSE), ]
    if (is.matrix(time_points)) data[[time_var]] <- c(t(time_points))
    data
}

extract_b <- function (object, id, n) {
    b <- data.matrix(ranef(object))
    mat <- matrix(0.0, n, ncol(b))
    colnames(mat) <- colnames(b)
    mat[id, ] <- b
    mat
}

extract_log_sigmas <- function (object) {
    out <- if (inherits(object, "lme")) {
        # we extract the log of sigma to be consistent with GLMMadaptive
        log(object$sigma)
    } else {
        object$phis
    }
    if (is.null(out)) out <- -20.0
    out
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

create_HC_X2 <- function (x, z, id) {
    check_tv <- function (x, id) {
        !all(sapply(split(x, id),
                    function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
    }

    cnams_x <- colnames(x)
    cnams_z <- colnames(z)
    if (!"(Intercept)" %in% cnams_x || !"(Intercept)" %in% cnams_z) {
        stop("cannot perform hierarchical centering in the absense of an ",
             "intercept term in both the fixed and random effects design ",
             "matrices.")
    }

    x_in_z <- which(cnams_x %in% cnams_z)
    x_notin_z <- which(!cnams_x %in% cnams_z)
    baseline <- x_notin_z[!apply(x[, x_notin_z, drop = FALSE], 2L, check_tv, id= id)]
    x_notin_z <- setdiff(x_notin_z, baseline)
    if (!length(baseline)) baseline <- as.integer(NA)
    if (!length(x_notin_z)) x_notin_z <- as.integer(NA)
    list(baseline = baseline, x_in_z = x_in_z, x_notin_z = x_notin_z,
         Xbase = x[!duplicated(id), baseline, drop = FALSE])
}

create_HC_X3 <- function(x, z, id, terms_FE, data, center = FALSE) {
    check_tv <- function (x, id) {
        !all(sapply(split(x, id),
                    function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
    }
    x <- scale(x, center = center, scale = FALSE)
    cnams_x <- colnames(x)
    cnams_z <- colnames(z)
    n_res <- ncol(z)
    X_HC <- vector("list", length = n_res)
    mat_HC <- matrix(0, nrow= n_res, ncol = ncol(x),
                     dimnames = list(cnams_z, cnams_x))
    mat_HC[cbind(which(cnams_z %in% cnams_x), which(cnams_x %in% cnams_z))] <- 1 # x_in_z
    # baseline (assumes every model has a random intercept)
    x_notin_z <- which(!cnams_x %in% cnams_z)
    ind <- !apply(x[, x_notin_z, drop = FALSE], 2L, check_tv, id = id)
    if(any(ind)) {
        baseline <- x_notin_z[ind]
        X_HC[[1]] <- x[!duplicated(id), baseline, drop = FALSE]
        mat_HC[cbind(1, baseline)] <- 2 # baseline
        }
    # remaining RE
    if (n_res > 1) {
        for (i in seq_len(n_res)[-1]) {
            xint_in_z <- union(grep(paste0(cnams_z[i], ":"), cnams_x, fixed= TRUE),
                               grep(paste0(":", cnams_z[i]), cnams_x, fixed= TRUE)) # interactions can be found as RE:var1, var1:RE, or var1:RE:var2
            if (!length(xint_in_z)) next
            data_temp <- data
            col_name <- colnames(data)[sapply(colnames(data), grepl, cnams_z[i], fixed= TRUE)]
            data_temp[, col_name][] <- 1
            x_temp <- scale(model.matrix.default(terms_FE, data = data_temp),
                            center = center, scale = FALSE)
            ind <- !apply(x_temp[, xint_in_z, drop = FALSE], 2L, check_tv,
                          id = id)
            if(any(ind)) {
                baseline_i <- xint_in_z[ind]
                X_HC[[i]] <- x_temp[!duplicated(id), baseline_i]
                mat_HC[cbind(i, baseline_i)] <- 3 # xint_in_z
            }
        }
    }
    x_in_z_base <- which(colSums(mat_HC > 0) == 1)
    x_notin_z <- which(colSums(mat_HC) == 0) # vars for MH
    if (!length(x_notin_z)) x_notin_z <- as.integer(NA)
    list(mat_HC = mat_HC, X_HC = X_HC, x_in_z_base = x_in_z_base,
         nfes_HC = length(x_in_z_base), z_in_x = which(rowSums(mat_HC == 1) == 1),
         x_in_z = which(colSums(mat_HC == 1) == 1), x_notin_z = x_notin_z,
         xbas_in_z = mat_HC[, x_in_z_base, drop = FALSE] > 1)
}

create_X_dot <- function(Xbase, nT, unq_idL, nres, nfes_HC, baseline, x_in_z_base, x_in_z) {

    n_outcomes <- length(nres) # number of outcomes
    n_res <- sum(nres) # total number of RE

    rows <- split(seq_len(n_res), rep(seq_along(nres), nres)) # all rows (id= 1)

    base_rows <- sapply(rows, head, 1) # rows for baseline (id= 1)
    base_cols <- mapply(function(xzb, b){ which(xzb %in% b)}, x_in_z_base, baseline) # cols for baseline

    RE_rows <- sapply(x_in_z, seq_along) # rows for RE (id= 1)
    RE_cols <- x_in_z # cols for RE

    M <- matrix(0, nrow= n_res*nT, ncol= sum(nfes_HC))

    for (j in seq_len(n_outcomes)) {

        ids <- unq_idL[[j]] # ids present in outcome-j
        ids_rows <- (ids-1) * n_res # 1st row for each id

        M[base_rows[j] + ids_rows, sum(nfes_HC[1:j-1]) + base_cols[[j]]] <- Xbase[[j]] # add baseline

        rows <- sum(nres[1:j-1]) + RE_rows[[j]] + rep(ids_rows, each= length(RE_rows[[j]]))
        cols <- rep(sum(nfes_HC[1:j-1]) + RE_cols[[j]], times= length(ids))
        M[cbind(rows, cols)] <- 1 # add 1 for each RE present in the FE
    }
    M
}

create_X_dot2 <- function (nT, nres, ind_FE_HC, x_in_z, x_in_z_base, unq_idL,
                           Xbase) {
    n_outcomes <- length (nres)
    ind_rows_subject <- rep(seq_len(nT), each = sum(nres))
    ind_rows_outcome <- rep(seq_len(n_outcomes), nres)
    ind_cols <- split(seq_along(ind_FE_HC),
                      rep(seq_len(n_outcomes), sapply(x_in_z_base, length)))

    M <- matrix(0.0, sum(nT * nres), length(ind_FE_HC))
    for (i in seq_len(nT)) {
        for (j in seq_len(n_outcomes)) {
            check <- i %in% unq_idL[[j]]
            if (check) {
                rows <- which(ind_rows_subject == i)[ind_rows_outcome == j]
                cols <- ind_cols[[j]][x_in_z[[j]]]
                M[cbind(rows[1:length(cols)], cols)] <- 1
                if (length(ind_cols[[j]]) > length(cols)) {
                    M[rows[1L], ind_cols[[j]][-x_in_z[[j]]]] <-
                        Xbase[[j]][as.character(i), ]
                }
            }
        }
    }
    M
}

create_X_dot3 <- function(nres, nfes_HC, z_in_x, x_in_z, X_HC, nT, unq_idL, xbas_in_z) {
    n_outcomes <- length(nres) # number of outcomes
    n_res <- sum(nres) # total number of RE
    M <- matrix(0, nrow = n_res * nT, ncol = sum(nfes_HC))
    for (j in seq_len(n_outcomes)) { # j-th outcome
        ids <- unq_idL[[j]] # ids present in outcome-j
        ids_rows <- (ids-1) * n_res # 1st row for each id
        rows <- sum(nres[1:j-1]) + z_in_x[[j]] + rep(ids_rows, each= length(z_in_x[[j]]))
        cols <- rep(sum(nfes_HC[1:j-1]) + x_in_z[[j]], times= length(ids))
        M[cbind(rows, cols)] <- 1 # add 1 for each z_in_x
        bas_cols <- xbas_in_z[[j]]
        for (k in z_in_x[[j]]) { # k-th RE in z_in_x # improve: remove this loop here
            if (!sum(bas_cols[k, ])) next
            M[sum(nres[1:j-1]) + k + ids_rows, sum(nfes_HC[1:j-1]) +
                  which(bas_cols[k, ])] <- X_HC[[j]][[k]]
        }
    }
    M
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
        unq_id_k <- unique(id_k)
        n <- length(unq_id_k)
        cov_postRE <- vector("list", n)
        names(cov_postRE) <- unq_id_k
        for (i in seq_len(n)) {
            Z_k_i <- Z_k[id_k == unq_id_k[i], , drop = FALSE]
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
    strata <- Data$strata
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
    n <- Data$n
    ###
    dataL <- data$dataL
    dataS <- data$dataS
    ###
    idVar <- model_info$var_names$idVar
    time_var <- model_info$var_names$time_var
    idT <- Data$idT
    terms_FE_noResp <- model_info$terms$terms_FE_noResp
    terms_RE <- model_info$terms$terms_RE
    terms_RE <- model_info$terms$terms_RE
    terms_Surv_noResp <- model_info$terms$terms_Surv_noResp
    dataL[[idVar]] <- dataL[[idVar]][drop = TRUE]
    ###
    functional_forms <- model_info$functional_forms
    FunForms_per_outcome <- model_info$FunForms_per_outcome
    collapsed_functional_forms <- model_info$collapsed_functional_forms
    ###
    Xbar <- Data$Xbar
    X_H <- Data$X_H; Z_H <- Data$Z_H; U_H <- Data$U_H; W0_H <- Data$W0_H; W_H <- Data$W_H
    X_h <- Data$X_h; Z_h <- Data$Z_h; U_h <- Data$U_h; W0_h <- Data$W0_h; W_h <- Data$W_h
    X_H2 <- Data$X_H2; Z_H2 <- Data$Z_H2; U_H2 <- Data$U_H2; W0_H2 <- Data$W0_H2; W_H2 <- Data$W_H2
    ###
    log_Pwk <- Data[["log_Pwk"]]
    log_Pwk2 <- Data[["log_Pwk2"]]
    ######################################################################################
    ######################################################################################
    times_long <- split(dataL[[time_var]], dataL[[idVar]])
    times_long <- times_long[sapply(times_long, length) > 0]
    dataS_init <- SurvData_HazardModel(times_long, dataS, Time_start,
                                         paste(idT, "_", strata), time_var)
    mf <- model.frame.default(terms_Surv_noResp, data = dataS_init)
    W_init <- construct_Wmat(terms_Surv_noResp, mf)
    if (!ncol(W_init)) {
        W_init <- cbind(W_init, rep(0, nrow(W_init)))
    }
    X_init <- design_matrices_functional_forms(times_long, terms_FE_noResp,
                                               dataL, time_var, idVar,
                                               collapsed_functional_forms, Xbar)
    Z_init <- design_matrices_functional_forms(times_long, terms_RE,
                                               dataL, time_var, idVar,
                                               collapsed_functional_forms)
    U_init <- lapply(functional_forms, construct_Umat, dataS = dataS_init)
    ##############
    fid <- dataL[[idVar]]
    fid <- factor(fid, levels = unique(fid))
    ind_multipl_events <-
        unlist(mapply2(rep.int,
                       x = lapply(split(fid, fid), function (x) as.numeric(unclass(x))),
                       times = tapply(idT, idT, length)), use.names = FALSE)
    id_init <- rep(list(dataL[[idVar]]), length.out = length(X_init))
    eta_init <- linpred_surv(X_init, betas, Z_init, b, id_init)
    eta_init[] <- mapply2(function (eta, ind) {
        if (length(ind) > length(eta)) eta[ind, , drop = FALSE] else eta
    }, eta_init, ind_multipl_events)
    Wlong_init <- create_Wlong(eta_init, FunForms_per_outcome, U_init)
    Wlong_init <- do.call("cbind", Wlong_init)
    ######################################################################################
    ######################################################################################
    start <- dataL[[time_var]]
    spl_Time <- split(Time_right, idT)
    stop <- unlist(mapply2(`c`, tapply(start, fid, tail, n = -1),
                           lapply(spl_Time, "[", 1L)), use.names = FALSE)
    create_event <- function (ni, delta) {
        if (ni == 1) delta else c(rep(0, ni - 1), delta)
    }
    event <- unlist(mapply2(create_event, ni = tapply(fid, fid, length),
                            delta = sapply(split(delta, idT), "[", 1L)), # <------ only one event, doesn't work with Comp Risks
                    use.names = FALSE)
    any_gammas <- !(ncol(W_init) == 1 && all(W_init[, 1] == 0))
    ind_multipl_events <-
        unlist(mapply2(rep.int,
                       x = lapply(split(fid, fid), unclass),
                       times = tapply(idT, idT, length)), use.names = FALSE)
    start <- start[ind_multipl_events]
    stop <- stop[ind_multipl_events]
    event <- event[ind_multipl_events]
    Wlong_init <- Wlong_init[ind_multipl_events, , drop = FALSE]
    WW <- if (any_gammas) cbind(W_init, Wlong_init) else Wlong_init
    ####
    fm <- coxph(Surv(start, stop, event) ~ WW)
    coefs <- coef(fm)
    gammas <- if (any_gammas) head(coefs, ncol(W_init)) else 0.0
    alphas <- tail(coefs, ncol(Wlong_init))
    alphas <- split(alphas, rep(seq_along(U_H), sapply(U_H, ncol)))
    V <- vcov(fm)
    if (any_gammas) {
        vcov_prop_gammas <- V[1:ncol(W_init), 1:ncol(W_init), drop = FALSE]
        vcov_prop_alphas <- V[-(1:ncol(W_init)), -(1:ncol(W_init)), drop = FALSE]
    } else {
        vcov_prop_gammas <- matrix(0.0, 1, 1)
        vcov_prop_alphas <- V
    }
    out <- list(gammas = gammas, alphas = alphas,
                vcov_prop_gammas = vcov_prop_gammas,
                vcov_prop_alphas = vcov_prop_alphas)
    ######################################################################################
    ######################################################################################
    id_H <- rep(list(rep(unclass(Data$idT), each = control$GK_k)), length(X_H))
    eta_H <- linpred_surv(X_H, betas, Z_H, b, id_H)
    Wlong_H <- create_Wlong(eta_H, FunForms_per_outcome, U_H)
    if (length(which_event)) {
        id_h <- rep(list(unclass(Data$idT)), length(X_h))
        eta_h <- linpred_surv(X_h, betas, Z_h, b, id_h)
        Wlong_h <- create_Wlong(eta_h, FunForms_per_outcome, U_h)
    } else {
        Wlong_h <- rep(list(matrix(0.0, length(Time_right), 1)), length(W_H))
    }
    if (length(which_interval)) {
        eta_H2 <- linpred_surv(X_H2, betas, Z_H, b, id_H)
        Wlong_H2 <- create_Wlong(eta_H2, FunForms_per_outcome, U_H2)
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
            H2 <- rowsum(exp(log_Pwk2 + lambda_H2), group = id_H[[1]], reorder = FALSE)
            log_Lik_surv[which_interval] <-
                log(exp(- H[which_interval]) -
                        exp(- (H2[which_interval] + H[which_interval])))
        }
        - sum(log_Lik_surv, na.rm = TRUE)
    }
    ncw <- ncol(W0_H)
    if(ncw > 20) {
        c(out, list(bs_gammas = rep(-0.1, ncw), vcov_prop_bs_gammas = diag(ncw)))
    } else {
        opt <- try({
            optim(rep(-0.1, ncw), log_dens_surv, method = "BFGS", hessian = TRUE,
                  control = list(parscale = rep(0.01, ncw)))
        }, silent = TRUE)
        V <- try(solve(opt$hessian), silent = TRUE)
        if (!inherits(opt, "try-error") && !inherits(V, "try-error")) {
            c(out, list(bs_gammas = opt$par, vcov_prop_bs_gammas = V))
        } else {
            c(out, list(bs_gammas = rep(-0.1, ncw),
                        vcov_prop_bs_gammas = diag(ncw)))
        }
    }
}

linpred_surv <- function (X, betas, Z, b, id) {
    out <- vector("list", length(X))
    for (i in seq_along(X)) {
        X_i <- X[[i]]
        Z_i <- Z[[i]]
        betas_i <- betas[[i]]
        b_i <- b[[i]]
        id_i <- id[[i]]
        out[[i]] <- matrix(0.0, nrow = nrow(X_i[[1]]), ncol = length(X_i))
        for (j in seq_along(X_i)) {
            X_ij <- X_i[[j]]
            Z_ij <- Z_i[[j]]
            out[[i]][, j] <- X_ij %*% betas_i + rowSums(Z_ij * b_i[id_i, ])
        }
    }
    out
}

FunForms_ind <- function (FunForms) {
    f <- function (l) rep(seq_along(l), sapply(l, length))
    lapply(FunForms, f)
}

create_Wlong <- function (eta, functional_forms, U) {
    Wlong <- vector("list", length(eta))
    for (i in seq_along(functional_forms)) {
        FF_i <- functional_forms[[i]]
        eta_i <- eta[[i]]
        U_i <- U[[i]]
        Wlong_i <- matrix(1.0, nrow(eta_i), max(unlist(FF_i)))
        for (j in seq_along(FF_i)) {
            ind <- FF_i[[j]]
            Wlong_i[, ind] <- Wlong_i[, ind] * eta_i[, j]
        }
        Wlong[[i]] <- U_i * Wlong_i
    }
    Wlong
}

create_Wlong2 <- function (eta, FunForms, U, ind) {
    Wlong <- vector("list", length(eta))
    for (i in seq_along(FunForms)) {
        FF_i <- FunForms[[i]]
        eta_i <- eta[[i]]
        U_i <- U[[i]]
        ind_i <- ind[[i]]
        Wlong_i <- matrix(1.0, nrow(eta_i), max(unlist(FF_i)))
        for (j in unique(ind_i)) {
            kk <- FF_i[ind_i == j]
            Wlong_i[, kk] <- Wlong_i[, kk] * eta_i[, j]
        }

        Wlong[[i]] <- U_i * Wlong_i
    }
    Wlong
}

Ptail <- function (x) {
    above <- mean(x >= 0)
    below <- mean(x < 0)
    2 * min(above, below)
}

effective_size <- function (x) {
    spectrum0.ar <- function(x) {
        d <- dim(x)
        nrx <- d[1L]
        ncx <- d[2L]
        v0 <- numeric(ncx)
        res <- as.matrix(lm.fit(cbind(1, seq_len(nrx)),
                                cbind(x, x))$residuals)
        for (i in seq_len(ncx)) {
            if (identical(all.equal(sd(res[, i]), 0), TRUE)) {
                v0[i] <- 0
            }
            else {
                ar.out <- ar(x[, i], aic = TRUE)
                v0[i] <- ar.out$var.pred / (1 - sum(ar.out$ar))^2
            }
        }
        v0
    }
    x <- as.matrix(x)
    spec <- spectrum0.ar(x)
    ifelse(spec == 0, 0, nrow(x) * apply(x, 2L, var) / spec)
}

std_err <- function (x) {
    x <- as.matrix(x)
    vars <- apply(x, 2L, var)
    ess <- effective_size(x)
    sqrt(vars / ess)
}

quantile2 <- function (x) quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)

cor2cov <- function (R, vars, sds = NULL) {
    p <- nrow(R)
    if (is.null(sds)) sds <- sqrt(vars)
    sds * R * rep(sds, each = p)
}

reconstr_D <- function (L, sds) {
    p <- length(sds)
    LL <- matrix(0.0, p, p)
    LL[upper.tri(LL)] <- L
    LL[1, 1] <- 1
    LL[cbind(2:p, 2:p)] <- sqrt(1 - colSums(LL^2)[-1L])
    out <- cor2cov(crossprod(LL), sds = sds)
    out[lower.tri(out, TRUE)]
}

lowertri2mat <- function (x, nams = NULL) {
    nx <- length(x)
    p <- round(0.5 * (sqrt(1 + 8 * nx) - 1))
    out <- matrix(0.0, p, p)
    out[lower.tri(out, TRUE)] <- x
    out[upper.tri(out)] <- t(out)[upper.tri(out)]
    out <- (out + t(out)) / 2
    if (!is.null(nams)) dimnames(out) <- list(nams, nams)
    out
}

lapply_nams <- function (X, FUN, ...) {
    out <- lapply(X, FUN, ...)
    names(out) <- X
    out
}

get_statistic <- function (s, stat) {
    out <- if (stat %in% c("Mean", "SD", "Time-series SE")) {
        s <- s$statistics
        if (is.matrix(s)) s[, stat] else s[stat]
    } else {
        s <- s$quantiles
        stat <- switch(stat, "Median" = "50%", "2.5CI" = "2.5%",
                       "97.5CI" = "97.5%")
        if (is.matrix(s)) s[, stat] else s[stat]
    }
}

center_fun <- function (M, means) {
    if (!all(M == 0)) as.matrix(M - rep(means, each = nrow(M))) else M
}

docall_cbind <- function (l) {
    if (is.list(l)) do.call("cbind", l) else l
}

printCall <- function (call) {
    d <- deparse(call)
    if (length(d) <= 3) {
        paste(d, sep = "\n", collapse = "\n")
    }
    else {
        d <- d[1:3]
        d[3] <- paste0(d[3], "...")
        paste(d, sep = "\n", collapse = "\n")
    }
}

mapply2 <- function (FUN, ..., MoreArgs = NULL, USE.NAMES = TRUE) {
    mapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = FALSE,
           USE.NAMES = USE.NAMES)
}

#help functions for ggplot mcmc diagnostics

# help function to extract mcmc lists
ggextractmcmc <- function(mcmc_list) {
    fun1 <- function(x) do.call(rbind, x)
    tmp <- lapply(mcmc_list, fun1)
    tmp2 <- lapply(mcmc_list, FUN = function(x) ncol(fun1(x)))
    list(do.call(cbind, tmp), do.call(c, tmp2))
}

# prepare data in a nice format to work with ggplot
# this can be exported in case a user wants to work with ggplot
# and use his/her own colors themes etc.
ggprepare <- function(object,
                      parm = c("all", "betas", "sigmas", "D", "bs_gammas",
                               "tau_bs_gammas", "gammas", "alphas")) {
    parm <- match.arg(parm)
    n_chains <- object$control$n_chains
    n_iter <- object$control$n_iter - object$control$n_burnin
    widedat_list <- ggextractmcmc(object$mcmc)
    widedat <- widedat_list[[1]]
    n_parms_each_fam <- widedat_list[[2]]
    n_parms <- ncol(widedat)
    parms <- colnames(widedat)
    parms <- make.unique(parms)
    parm_fam <- names(object$mcmc)
    reps <- n_parms_each_fam * (n_iter * n_chains)
    parm_fam <- rep(parm_fam, times = reps)
    parm_fam <- gsub('[[:digit:]]+', '', parm_fam)
    ggdata <- expand.grid('iteration' = 1:n_iter,
                          'chain' = 1:n_chains,
                          'parm' = parms)
    ggdata$value <- as.vector(widedat)
    ggdata$parm_fam <- parm_fam
    ggdata <- ggdata[, c('iteration', 'chain', 'parm_fam', 'parm', 'value')]
    ggdata$chain <- factor(ggdata$chain)
    ggdata$parm_fam <- factor(ggdata$parm_fam, levels = unique(parm_fam))
    ggdata$parm <- factor(ggdata$parm, levels = unique(ggdata$parm))
    if (parm == "all") {
        ggdata
    } else {
        ggdata[ggdata$parm_fam %in% parm, ]
    }
}

# fancy color themes
# This is not a function and is better to not be exported but be
# an object used only internally
ggcolthemes <- list(
    'standard' = c("1" = '#363636', "2" = '#f25f5c', "3" = '#247ba0'),
    'catalog' = c("1" = '#cc2a36', "2" = '#edc951', "3" = '#00a0b0'),
    'metro' = c("1" = '#d11141', "2" = '#00aedb', "3" = '#ffc425'),
    'pastel' = c("1" = '#a8e6cf', "2" = '#ffd3b6', "3" = '#ff8b94'),
    'beach' = c("1" = '#ff6f69', "2" = '#ffcc5c', "3" = '#88d8b0'),
    'moonlight' = c("1" = '#3da4ab', "2" = '#f6cd61', "3" = '#fe8a71'),
    'goo' = c("1" = '#008744', "2" = '#0057e7', "3" = '#d62d20'),
    'sunset' = c("1" = '#f67e7d', "2" = '#843b62', "3" = '#0b032d')
)

tv <- function (x, knots = NULL, ord = 2L) {
    if (is.null(knots)) {
        knots <- .knots_base_hazard
    }
    out <- splines::splineDesign(knots, x, ord = ord, outer.ok = TRUE)
    attr(out, "knots") <- knots
    attr(out, "ord") <- ord
    attr(out, "class") <- c("tve", "basis", "matrix")
    out
}

makepredictcall.tv <- function (var, call) {
    if (as.character(call)[1L] != "tve")
        return(call)
    at <- attributes(var)[c("knots", "ord")]
    x <- call[1L:2L]
    x[names(at)] <- at
    x
}

create_W0 <- function (times, knots, ord, strata) {
    W0 <- splineDesign(knots, times, ord = ord, outer.ok = TRUE)
    n_strata <- length(unique(strata))
    ncW0 <- ncol(W0)
    ind_cols <- matrix(seq_len(ncW0 * n_strata), ncW0)
    out <- matrix(0.0, nrow(W0), ncW0 * n_strata)
    for (i in seq_len(n_strata)) {
        row_inds <- strata == i
        col_inds <- ind_cols[, i]
        out[row_inds, col_inds] <- W0[row_inds, ]
    }
    out
}

construct_Umat <- function (fForms, dataS) {
    tt <- terms(fForms)
    m <- model.matrix(tt, model.frame(tt, data = dataS))
    cnams <- colnames(m)
    ind_value <- grep("value(", cnams, fixed = TRUE)
    ind_slope <- grep("slope(", cnams, fixed = TRUE)
    ind_area <- grep("area(", cnams, fixed = TRUE)
    ind <- unique(c(ind_value, ind_slope, ind_area))
    m[, cnams %in% cnams[ind], drop = FALSE]
}

construct_Wmat <- function (Terms, mf) {
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
            if (any(attr(Terms, "order")[attr(Terms, "factors")[i, ] > 0] > 1))
                hasinteractions <- TRUE
        }
        if (!hasinteractions)
            dropterms <- stemp$terms
    }
    if (length(dropterms)) {
        Terms2 <- Terms[-dropterms]
        X <- model.matrix(Terms2, mf)
        temp <- attr(X, "assign")
        shift <- sort(dropterms)
        temp <- temp + 1 * (shift[1] <= temp)
        if (length(shift) == 2)
            temp + 1 * (shift[2] <= temp)
        attr(X, "assign") <- temp
    } else X <- model.matrix(Terms, mf)
    Xatt <- attributes(X)
    adrop <- if (hasinteractions) c(0, untangle.specials(Terms, "strata")$terms) else 0
    xdrop <- Xatt$assign %in% adrop
    X <- X[, !xdrop, drop = FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
    attr(X, "contrasts") <- Xatt$contrasts
    X
}

nearPD <- function (M, eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08,
                    maxits = 100) {
    if (!is.numeric(M))
        stop("Input matrix 'M' must be numeric.")
    if (length(M) == 1)
        return(abs(M))
    if (is.matrix(M) && !identical(M, t(M)))
        stop("Input matrix M must be square and symmetric.\n")
    inorm <- function (x) max(rowSums(abs(x)))
    n <- ncol(M)
    U <- matrix(0.0, n, n)
    X <- M
    iter <- 0
    converged <- FALSE
    while (iter < maxits && !converged) {
        Y <- X
        T <- Y - U
        e <- eigen(Y, symmetric = TRUE)
        Q <- e$vectors
        d <- e$values
        D <- if (length(d) > 1) diag(d) else as.matrix(d)
        p <- (d > eig.tol * d[1])
        QQ <- Q[, p, drop = FALSE]
        X <- QQ %*% D[p, p, drop = FALSE] %*% t(QQ)
        U <- X - T
        X <- (X + t(X)) / 2
        conv <- inorm(Y - X) / inorm(Y)
        iter <- iter + 1
        converged <- conv <= conv.tol
    }
    X <- (X + t(X)) / 2
    e <- eigen(X, symmetric = TRUE)
    d <- e$values
    Eps <- posd.tol * abs(d[1])
    if (d[n] < Eps) {
        d[d < Eps] <- Eps
        Q <- e$vectors
        o.diag <- diag(X)
        X <- Q %*% (d * t(Q))
        D <- sqrt(pmax(Eps, o.diag) / diag(X))
        X[] <- D * X * rep(D, each = n)
    }
    (X + t(X)) / 2
}

fit_stats <- function (lL, lL_mean_parms) {
    D_bar <- - 2.0 * mean(rowSums(lL, na.rm = TRUE))
    D_hat <- - 2.0 * sum(lL_mean_parms, na.rm = TRUE)
    pD <- D_bar - D_hat
    CPO <- 1 / colMeans(exp(-lL), na.rm = TRUE)
    CC <- log(nrow(lL))
    LPML <- sum(- colLogSumExps(-lL, na.rm = TRUE) + CC)
    LPPD <- - 2.0 * sum(colLogSumExps(lL, na.rm = TRUE) - CC)
    pWAIC2 <- 2.0 * sum(colVars(lL, na.rm = TRUE))
    list(DIC = pD + D_bar, pD = pD, LPML = LPML,
         CPO = CPO, WAIC = LPPD + pWAIC2)
}

get_betas_nHC <- function (v, ind) {
    if (any(is.na(ind))) {
        if (is.matrix(v)) matrix(1.0) else 0.0
    } else {
        if (is.matrix(v)) v[ind, ind, drop = FALSE] else v[ind]
    }
}

weak_informative_Tau <- function (model, Xbar) {
    V <- vcov_center(vcov2(model), Xbar)
    diags <- 150.0 * diag(V)
    diag(1 / diags, nrow(V), ncol(V))
}

weak_informative_Tau2 <- function (y, X, is_gaussian) {
    s_y <- if (is_gaussian) sd(y) else 1.0
    s_x <- apply(X, 2L, sd)
    if (colnames(X)[1L] == "(Intercept)") s_x[1L] <- 1
    diag(s_x / s_y, length(s_x)) / 30
}

weak_informative_mean <- function (y, X, is_gaussian) {
    out <- rep(0.0, ncol(X))
    #out[1L] <- if (is_gaussian) mean(y) else 0.0
    out
}

vcov_center <- function (vcov, Xbar) {
    if (!all(abs(Xbar) < sqrt(.Machine$double.eps))) {
        Xbar[1L] <- 1
        var_temp <- crossprod(Xbar, vcov) %*% Xbar
        vcov[1L, -1L] <- vcov[-1L, 1L] <- c(vcov %*% Xbar)[-1L]
        vcov[1L, 1L] <- var_temp
    }
    vcov
}

center_X <- function (x, ind) {
    if (length(ind) > 1 || !is.na(ind)) x[-ind] <- 0.0 else x <- x * 0
    x
}

jitter2 <- function (x, factor = 2) {
    if (is.list(x)) {
        x[] <- lapply(x, jitter, factor = factor)
    } else {
        jitter(x, factor = factor)
    }
}



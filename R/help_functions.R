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
    sapply(c("value", "slope", "area"), grep, x = colnames(M), fixed = TRUE,
           simplify = FALSE)
}

expand_Dexps <- function (Form, respVar) {
    tlabs <- attr(terms(Form), "term.labels")
    dexps_ind <- grep("Dexp", tlabs)
    dexps <- tlabs[dexps_ind]
    dexps <- gsub("slope", "value", dexps)
    p1 <- paste0(respVar, "))")
    p2 <- paste0(respVar, ")):slope(", respVar, ")")
    dexps <- gsub(p1, p2, dexps, fixed = TRUE)
    tlabs[dexps_ind] <- dexps
    reformulate(tlabs)
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

extract_attributes <- function (form, data) {
    mf <- model.frame(terms(form), data = data)
    eps <- lapply(mf, function (v, name) attr(v, name), name = "eps")
    direction <- lapply(mf, function (v, name) attr(v, name), name = "direction")
    list(eps = eps[!sapply(eps, is.null)],
         direction = direction[!sapply(direction, is.null)])
}

design_matrices_functional_forms <-
    function (time, terms, data, timeVar, idVar, idT, Fun_Forms, Xbar = NULL,
              eps, direction) {
        data[] <- lapply(data, function (x) locf(locf(x), fromLast = TRUE))
        desgn_matr <- function (time, terms, Xbar) {
            D <- LongData_HazardModel(time, data, data[[timeVar]],
                                      data[[idVar]], idT, timeVar)
            mf <- lapply(terms, model.frame.default, data = D)
            X <- mapply2(model.matrix.default, terms, mf)
            if (!is.null(Xbar))
                X <- mapply2(function (m, mu) m - rep(mu, each = nrow(m)), X, Xbar)
            X
        }
        degn_matr_slp <- function (time, terms, Xbar, eps, direction) {
            K <- length(terms)
            out <- vector("list", K)
            for (i in seq_len(K)) {
                direction_i <- if (length(direction[[i]])) direction[[i]][[1L]] else "both"
                eps_i <- if (length(eps[[i]])) eps[[i]][[1L]] else 0.001
                if (direction_i == "both") {
                    t1 <- time + eps_i
                    t2 <- pmax(time - eps_i, 0)
                    e <- 2 * eps_i
                } else {
                    t1 <- time
                    t2 <- pmax(time - eps_i, 0)
                    e <- eps_i
                }
                terms_i <- terms[[i]]
                D1 <- LongData_HazardModel(t1, data, data[[timeVar]],
                                           data[[idVar]], idT, timeVar)
                mf1 <- model.frame.default(terms_i, data = D1)
                X1 <- model.matrix.default(terms_i, mf1)
                D2 <- LongData_HazardModel(t2, data, data[[timeVar]],
                                           data[[idVar]], idT, timeVar)
                mf2 <- model.frame.default(terms_i, data = D2)
                X2 <- model.matrix.default(terms_i, mf2)
                if (!is.null(Xbar)) {
                    X1 <- X1 - rep(Xbar[[i]], each = nrow(X1))
                    X2 <- X2 - rep(Xbar[[i]], each = nrow(X2))
                }
                out[[i]] <- (X1 - X2) / e
            }
            out
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
                # we divide with x to obtain the area up to time t, divided by t
                # to account for the length of the interval
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
                    "slope" = degn_matr_slp(time, terms, Xbar, eps, direction),
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

value <- area <- function (x) rep(1, NROW(x))
vexpit <- Dexpit <- vexp <- Dexp <- function (x) rep(1, NROW(x))
vsqrt <- vlog <- vlog2 <- vlog10 <- function (x) rep(1, NROW(x))
slope <- function (x, eps = 0.001, direction = "both") {
    out <- rep(1, NROW(x))
    temp <- list(eps = eps, direction = direction)
    attributes(out) <- c(attributes(out), temp)
    out
}

create_HC_X <- function(x, z, id, terms, data) {
    check_tv <- function (x, id) {
        !all(sapply(split(x, id),
                    function (z) all(abs(z - z[1L]) < .Machine$double.eps^0.5)))
    }
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
            # interactions can be found as RE:var1, var1:RE, or var1:RE:var2
            xint_in_z <- union(grep(paste0(cnams_z[i], ":"), cnams_x,
                                    fixed = TRUE),
                               grep(paste0(":", cnams_z[i]), cnams_x,
                                    fixed = TRUE))
            xint_in_z <- sort(xint_in_z)
            if (!length(xint_in_z)) next
            data_temp <- data
            col_name <- colnames(data)[sapply(colnames(data), grepl, cnams_z[i],
                                              fixed = TRUE)]
            data_temp[, col_name][] <- 1
            x_temp <- model.matrix.default(terms, data = data_temp)
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

create_X_dot <- function(nres, nfes_HC, z_in_x, x_in_z, X_HC, nT, unq_idL, xbas_in_z) {
    n_outcomes <- length(nres) # number of outcomes
    n_res <- sum(nres) # total number of RE
    M <- matrix(0, nrow = n_res * nT, ncol = sum(nfes_HC))
    for (j in seq_len(n_outcomes)) { # j-th outcome
        ids <- unq_idL[[j]] # ids present in outcome-j
        ids_rows <- (ids-1) * n_res # 1st row for each id
        rows1 <- sum(nres[1:j-1]) + z_in_x[[j]] + rep(ids_rows,
                                                      each = length(z_in_x[[j]]))
        cols1 <- sum(nfes_HC[1:j-1]) + match(names(x_in_z[[j]]),
                                             colnames(xbas_in_z[[j]]))
        cols1 <- rep(cols1, times = length(ids))
        M[cbind(rows1, cols1)] <- 1 # add 1 for each z_in_x
        bas_cols <- xbas_in_z[[j]]
        for (k in seq_len(nres[j])) { # k-th RE
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

extractFuns_FunForms <- function (Form, data) {
    tr <- terms(Form)
    mF <- model.frame(tr, data = data)
    M <- model.matrix(tr, mF)
    cnams <- colnames(M)
    possible_forms <- c("value(", "slope(", "area(")
    ind <- unlist(lapply(possible_forms, grep, x = cnams, fixed = TRUE))
    M <- M[1, cnams %in% cnams[unique(ind)], drop = FALSE]
    FForms <- sapply(c("value", "slope", "area"), grep, x = colnames(M),
                     fixed = TRUE, simplify = FALSE)
    FForms <- FForms[sapply(FForms, length) > 0]
    get_fun <- function (FForm, nam) {
        cnams <- colnames(M)[FForm]
        out <- rep("identity", length(cnams))
        f <- function (fun_nam) {
            grep(paste0(fun_nam, "(", nam), cnams, fixed = TRUE)
        }
        out[f("expit")] <- "expit"
        out[f("dexpit")] <- "dexpit"
        out[f("exp")] <- "exp"
        out[f("dexp")] <- "dexp"
        out[f("log")] <- "log"
        out[f("log2")] <- "log2"
        out[f("log10")] <- "log10"
        out[f("sqrt")] <- "sqrt"
        out[f("square")] <- "square"
        out[1L] # <- to change: if the same term different functions
    }
    mapply(get_fun, FForms, names(FForms))
}

transf_eta <- function (eta, fun_nams) {
    for (j in seq_along(fun_nams)) {
        if (fun_nams[j] == "expit") {
            eta[, j] <- plogis(eta[, j])
        } else if (fun_nams[j] == "dexpit") {
            eta[, j] <- plogis(eta[, j]) * plogis(eta[, j], lower.tail = FALSE)
        } else if (fun_nams[j] == "exp") {
            eta[, j] <- exp(eta[, j])
        } else if (fun_nams[j] == "log") {
            eta[, j] <- log(eta[, j])
        } else if (fun_nams[j] == "sqrt") {
            eta[, j] <- sqrt(eta[, j])
        } else if (fun_nams[j] == "square") {
            eta[, j] <- eta[, j] * eta[, j]
        }
    }
    eta
}

create_Wlong <- function (eta, functional_forms, U, Funs_FunsForms) {
    Wlong <- vector("list", length(eta))
    for (i in seq_along(functional_forms)) {
        FF_i <- functional_forms[[i]]
        eta_i <- transf_eta(eta[[i]], Funs_FunsForms[[i]])
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

center_fun <- function (M, means, sds = NULL) {
    if (!all(M == 0)) {
        if (is.null(sds)) {
            scale(x = M, center = means, scale = FALSE)
        } else {
            scale(x = M, center = means, scale = sds)
        }
    } else M
}

docall_cbind <- function (l) {
    if (is.list(l)) do.call("cbind", l) else l
}

docall_rbind <- function (l) {
    if (is.list(l)) do.call("rbind", l) else l
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
    mcmc_list <- mcmc_list[!names(mcmc_list) %in% 'b']
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
    n_iter <- (object$control$n_iter - object$control$n_burnin) / object$control$n_thin
    widedat_list <- ggextractmcmc(object$mcmc)
    widedat <- widedat_list[[1]]
    n_parms_each_fam <- widedat_list[[2]]
    n_parms <- ncol(widedat)
    parms <- colnames(widedat)
    parms <- make.unique(parms)
    parm_fam <- names(object$mcmc[!names(object$mcmc) %in% 'b'])
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
    W0 <- lapply(knots, splineDesign, x = times, ord = ord, outer.ok = TRUE)
    n_strata <- length(unique(strata))
    ncW0 <- ncol(W0[[1L]])
    ind_cols <- matrix(seq_len(ncW0 * n_strata), ncW0)
    out <- matrix(0.0, nrow(W0[[1L]]), ncW0 * n_strata)
    for (i in seq_len(n_strata)) {
        row_inds <- strata == i
        col_inds <- ind_cols[, i]
        out[row_inds, col_inds] <- W0[[i]][row_inds, ]
    }
    out
}

construct_Umat <- function (fForms, dataS) {
    #expit <- exp <- log <- log2 <- log10 <- sqrt <- function (x) rep(1, NROW(x))
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
    diags <- pmin(100.0 * diag(V), 10)
    #diags <- rep(100, length(diags))
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

SurvData_HazardModel <- function (time_points, data, times, ids, time_var) {
    unq_ids <- unique(ids)
    fids <- factor(ids, levels = unq_ids)
    if (!is.matrix(time_points)) {
        time_points <- as.matrix(time_points)
    }
    tt <- split(time_points, row(time_points))
    spl_times <- split(times, fids)[unclass(fids)]
    rownams_id <- split(row.names(data), fids)[unclass(fids)]
    ind <- mapply2(findInterval, tt, spl_times)
    ind[] <- lapply(ind, function (x) {x[x < 1] <- 1; x})
    ind <- mapply2(`[`, rownams_id, ind)
    data <- data[unlist(ind, use.names = FALSE), ]
    data[[time_var]] <- c(t(time_points))
    row.names(data) <- seq_len(nrow(data))
    data
}

LongData_HazardModel <- function (time_points, data, times, ids, idT, time_var) {
    unq_ids <- unique(ids)
    fids <- factor(ids, levels = unq_ids)
    idT. <- match(idT, unique(idT))
    tt <- if (is.list(time_points)) {
        time_points
    } else {
        if (!is.matrix(time_points)) {
            time_points <- as.matrix(time_points)
        }
        split(time_points, row(time_points))
    }
    spl_times <- split(times, fids)
    rownams_id <- split(row.names(data), fids)
    if (length(spl_times) < length(tt)) {
        spl_times <- spl_times[idT.]
        rownams_id <- rownams_id[idT.]
    }
    ind <- mapply2(findInterval, tt, spl_times)
    ind[] <- lapply(ind, function (x) {x[x < 1] <- 1; x})
    ind <- mapply2(`[`, rownams_id, ind)
    data <- data[unlist(ind, use.names = FALSE), ]
    data[[time_var]] <- if (is.matrix(time_points)) {
        c(t(time_points))
    } else {
        unlist(time_points, use.names = FALSE)
    }
    row.names(data) <- seq_len(nrow(data))
    data
}

ms_setup <- function (data, timevars, statusvars, transitionmat, id, covs = NULL) {
    # setup times matrix with NAs
    # First row is NA as this is starting state
    timesmat <- matrix(NA, nrow(data), length(timevars))
    timecols_data <- which(colnames(data) %in% timevars[!is.na(timevars)])
    timesmat[, -which(is.na(timevars))] <- as.matrix(data[, timecols_data])
    # setup status matrix with NAs
    # First row is NA as this is starting state
    statusmat <- matrix(NA, nrow(data), length(statusvars))
    statuscols_data <- which(colnames(data) %in% statusvars[!is.na(statusvars)])
    statusmat[, -which(is.na(statusvars))] <- as.matrix(data[, statuscols_data])
    # ensure convert to matrices
    timesmat <- as.matrix(timesmat)
    statusmat <- as.matrix(statusmat)
    # check dimesnions are the same
    if (any(dim(timesmat) != dim(statusmat)))
        stop("Dimensions of \"time\" and \"status\" data should be equal")
    # components
    # number of unique subjects
    n_subj <- nrow(timesmat)
    # number of states
    n_states <- dim(transitionmat)[1]
    # set start state to 1 and start time to 0 for all subjects
    # ATTENTION: this needs to be adjusted to more flexible to allow subjects starting at different states
    # this could be achieved by a requesting a separate argument (vector with starting state)
    starting_state <- rep(1, n_subj)
    starting_time <- rep(0, n_subj)
    idnam <- id
    id <- data[[id]]
    order_id <- order(id)
    out <- ms_prepdat(timesmat = timesmat, statusmat = statusmat, id = id,
                   starting_time = starting_time, starting_state = starting_state,
                   transitionmat = transitionmat,
                   original_states = (1:nrow(transitionmat)), longmat = NULL)
    out <- as.data.frame(out)
    names(out) <- c(idnam, "from_state", "to_state", "transition",
                    "Tstart", "Tstop", "status")
    out$time <- out$Tstop - out$Tstart
    out <- out[, c(1:6, 8, 7)]
    ord <- order(out[, 1], out[, 5], out[, 2], out[, 3])
    out <- out[ord, ]
    row.names(out) <- 1:nrow(out)
    # Covariates
    if (!is.null(covs)) {
        n_covs <- length(covs)
        cov_cols <- match(covs, names(data))
        cov_names <- covs
        covs <- data[, cov_cols]
        if (!is.factor(out[, 1]))
            out[, 1] <- factor(out[, 1])
        n_per_subject <- tapply(out[, 1], out[, 1], length)
        if (n_covs > 1)
            covs <- covs[order_id, , drop = FALSE]
        if (n_covs == 1) {
            longcovs <- rep(covs, n_per_subject)
            longcovs <- longcovs[ord]
            longcovs <- as.data.frame(longcovs)
            names(longcovs) <- cov_names
        } else {
            longcovs <- lapply(1:n_covs, function(i) rep(covs[, i], n_per_subject))
            longcovs <- as.data.frame(longcovs)
            names(longcovs) <- cov_names
        }
        out <- cbind(out, longcovs)
    }
    # add attributes maybe
    # add specific class maybe
    # need to add functionality for covariates (e.g. like keep in mstate)
    return(out)
}

ms_prepdat <- function (timesmat, statusmat, id, starting_time, starting_state, transitionmat,
                        original_states, longmat) {
    if (is.null(nrow(timesmat)))
        return(longmat)
    if (nrow(timesmat) == 0)
        return(longmat)
    from_states <- apply(!is.na(transitionmat), 2, sum)
    to_states <- apply(!is.na(transitionmat), 1, sum)
    absorbing_states <- which(to_states == 0)
    starts <- which(from_states == 0)
    new_states <- starting_state
    new_times <- starting_time
    rmv <- NULL
    for (i in 1:starts) {
        subjects <- which(starting_state == starts)
        n_start <- length(subjects)
        to_states_2 <- which(!is.na(transitionmat[starts, ]))
        trans_states <- transitionmat[starts, to_states_2]
        n_trans_states <- length(to_states_2)
        if (all(n_start > 0, n_trans_states > 0)) {
            Tstart <- starting_time[subjects]
            Tstop <- timesmat[subjects, to_states_2, drop = FALSE]
            Tstop[Tstop <= Tstart] <- Inf
            state_status <- statusmat[subjects, to_states_2, drop = FALSE]
            mintime <- apply(Tstop, 1, min)
            hlp <- Tstop * 1 / state_status
            hlp[Tstop == 0 & state_status == 0] <- Inf
            next_time <- apply(hlp, 1, min)
            censored <- which(is.infinite(apply(hlp, 1, min)))
            wh <- which(mintime < next_time)
            whminc <- setdiff(wh, censored)
            if (length(whminc) > 0) {
                whsubjs <- id[subjects[whminc]]
                whsubjs <- paste(whsubjs, collapse = " ")
                warning("Subjects ", whsubjs, " Have smaller transition time with status = 0, larger transition time with status = 1,
                from starting state ", original_states[starts])
            }
            next_time[censored] <- mintime[censored]
            if (ncol(hlp) > 1) {
                hlpsrt <- t(apply(hlp, 1, sort))
                warn1 <- which(hlpsrt[, 1] - hlpsrt[, 2] == 0)
                if (length(warn1) > 0) {
                    isw <- id[subjects[warn1]]
                    isw <- paste(isw, collapse = " ")
                    hsw <- hlpsrt[warn1, 1]
                    hsw <- paste(hsw, collapse = " ")
                    warning("simultaneous transitions possible for subjects ", isw, " at times ", hsw,
                            " -> Smallest receiving state will be used")
                }
            }
            if (length(censored) > 0) {
                next_state <- apply(hlp[-censored, , drop = FALSE],
                                    1, which.min)
                absorbed <- (1:n_start)[-censored][which(to_states_2[next_state] %in% absorbing_states)]
            } else {
                next_state <- apply(hlp, 1, which.min)
                absorbed <- (1:n_start)[which(to_states_2[next_state] %in% absorbing_states)]
            }
            states_matrix <- matrix(0, n_start, n_trans_states)
            if (length(censored) > 0) {
                states_matrix_min <- states_matrix[-censored, , drop = FALSE]
            } else {
                states_matrix_min <- states_matrix
            }
            if (nrow(states_matrix_min) > 0)
                states_matrix_min <- t(sapply(1:nrow(states_matrix_min), function(i) {
                    x <- states_matrix_min[i, ]
                    x[next_state[i]] <- 1
                    return(x)
                }))
            if (length(censored) > 0) {
                states_matrix[-censored, ] <- states_matrix_min
            } else {
                states_matrix <- states_matrix_min
            }
            mm <- matrix(c(rep(id[subjects], rep(n_trans_states, n_start)),
                           rep(original_states[starts], n_trans_states * n_start),
                           rep(original_states[to_states_2], n_start),
                           rep(trans_states, n_start), rep(Tstart, rep(n_trans_states, n_start)),
                           rep(next_time, rep(n_trans_states, n_start)), as.vector(t(states_matrix))),
                         n_trans_states * n_start, 7)
            longmat <- rbind(longmat, mm)
            rmv <- c(rmv, subjects[c(censored, absorbed)])
            if (length(censored) > 0) {
                new_states[subjects[-censored]] <- to_states_2[next_state]
            } else {
                new_states[subjects] <- to_states_2[next_state]
            }
            if (length(censored) > 0)  {
                new_times[subjects[-censored]] <- next_time[-censored]
            } else {
                new_times[subjects] <- next_time
            }
        }
    }
    if (length(rmv) > 0) {
        timesmat <- timesmat[-rmv, ]
        statusmat <- statusmat[-rmv, ]
        new_times <- new_times[-rmv]
        new_states <- new_states[-rmv]
        id <- id[-rmv]
    }
    n_states <- nrow(transitionmat)
    idx <- rep(1, n_states)
    idx[starts] <- 0
    idx <- cumsum(idx)
    new_states <- idx[new_states]
    Recall(timesmat = timesmat[, -starts], statusmat = statusmat[, -starts],
           id = id, starting_time = new_times, starting_state = new_states,
           transitionmat = transitionmat[-starts, -starts], original_states = original_states[-starts],
           longmat = longmat)
}

extract_mcmc_as_inits <- function(x_mcmc, i, nams_x, nams_special, ind_RE, dim_D, has_sigmas) {
    c(lapply(x_mcmc[!nams_x %in% nams_special], function(x, i) x[[i]][nrow(x[[i]]), , drop = TRUE], i = i),
      'log_sigmas' = unname(lapply(x_mcmc[nams_x %in% 'sigmas'], function(x, i, has_sigmas) ifelse(has_sigmas, log(x[[i]][nrow(x[[i]]), , drop = TRUE]), -20), i = i, has_sigmas = has_sigmas)),
      list(betas = unname(lapply(x_mcmc[nams_x[grep('betas', nams_x)]], function(x, i) x[[i]][nrow(x[[i]]), , drop = TRUE], i = i))),
      list(b = mapply(function(x, y, i) x[[i]][,y , dim(x[[i]])[length(dim(x[[i]]))], drop = TRUE],
                      x_mcmc[nams_x %in% 'b'], ind_RE,
                      i = i, SIMPLIFY = FALSE, USE.NAMES = FALSE)),
      lapply(x_mcmc[nams_x %in% 'D'], function(x, i, dim_D) {
          D <- matrix(NA, ncol = dim_D, nrow = dim_D)
          D_vec <- x[[i]][nrow(x[[i]]), , drop = TRUE]
          D[upper.tri(D, diag = TRUE)] <- D_vec[order(names(D_vec))]
          D[lower.tri(D, diag = FALSE)] <- D[upper.tri(D, diag = FALSE)]
          D
      },
      i = i, dim_D = dim_D)
    )
}

extract_last_iterations <- function (x) {
    x_mcmc <- x$mcmc
    n_chains <- x$control$n_chains
    nams_x <- names(x_mcmc)
    nams_special <- c(nams_x[grep('betas', nams_x)], 'b', 'D', 'sigmas')
    ind_RE <- x$model_data$ind_RE
    dim_D <- ncol(x$statistics$Mean$b)
    last_iter_x <- vector('list', n_chains)
    has_sigmas <- as.integer(x$initial_values$log_sigmas > -20)
    for (i in 1:n_chains) {
        last_iter_x[[i]] <-
            extract_mcmc_as_inits(x_mcmc, i = i, nams_x = nams_x,
                                  nams_special = nams_special, ind_RE = ind_RE,
                                  dim_D, has_sigmas = has_sigmas)
    }
    last_iter_x
}

create_sigma_list <- function (sigmas, ss_sigmas, idL) {
    n <- length(sigmas)
    out <- vector("list", n)
    for (i in 1:n) {
        sigmas_i <- sigmas[[i]]
        id_i <- idL[[i]]
        if (ss_sigmas[i]) {
            out[[i]] <- sigmas_i[id_i]
        } else {
            out[[i]] <- rep(sigmas_i, length(id_i))
        }
    }
    out
}

lng_unq <- function (x) length(unique(x))


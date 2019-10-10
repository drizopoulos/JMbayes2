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

right_rows <- function (data, times, ids, time_points) {
    # data = dataL
    # times = dataL[[timeVar]]
    # ids = dataL[[idVar]]
    # time_points = Time


    unq_ids <- unique(ids)
    fids <- factor(ids, levels = unq_ids)
    if (!is.matrix(time_points)) {
        time_points <- as.matrix(time_points)
    }
    if (nrow(time_points) != length(unq_ids)) {
        stop("the length of unique 'ids' does not match the number of rows ",
             "of 'time_points'.")
    }
    time_points <- split(time_points, row(time_points))
    ind <- mapply(findInterval, time_points, split(times, fids))
    ind[ind < 1] <- 1
    if (!is.matrix(ind)) {
        ind <- rbind(ind)
    }
    rownams_id <- split(row.names(data), fids)
    ind <- mapply(`[`, rownams_id, split(ind, col(ind)))
    data[c(ind), ]
}

last_rows <- function (data, ids) {
    fidVar <- factor(ids, levels = unique(ids))
    data[tapply(row.names(data), fidVar, tail, n = 1L), ]
}

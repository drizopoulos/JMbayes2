if (FALSE) {
    library("JMbayes2")
    pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
    CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
    fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
    fm2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)
    fm3 <- mixed_model(ascites ~ year + sex, data = pbc2,
                       random = ~ year | id, family = binomial())

    jointFit1 <- jm(CoxFit, list(fm1, fm2, fm3), time_var = "year")

    source("./R/help_functions.R")
}


object <- jointFit1
newdata <- pbc2[pbc2$id %in% c(2, 3, 15), ]
newdata$id <- factor(newdata$id)
newdata2 <- newdata[newdata$year > 1, ]
newdata <- newdata[newdata$year < 1, ]
newdata$status2 <- 0
newdata$years <- with(newdata, ave(year, id, FUN = function (x) max(x, na.rm = T)))
newdata$prothrombin[c(3, 5, 8)] <- NA
process <- c("Event")
level <- 1L
CI_level <- 0.95
pred_type <- "response"
return_newdata <- TRUE
n_samples <- 200L
n_mcmc <- 25L
cores <- 3L
seed <- 1L

#############################################################
#############################################################

components_newdata <- get_components_newdata(object, newdata)

predict_Long <- function (components_newdata, newdata2, level, pred_type,
                          CI_level) {
    # Predictions for newdata
    betas <- components_newdata$mcmc[["betas"]]
    b_mat <- components_newdata$mcmc[["b"]]
    ind_RE <- components_newdata$ind_RE
    links <- components_newdata$links
    K <- length(ind_RE)
    M <- dim(b_mat)[3L]
    out <- lapply(components_newdata$X, function (x) matrix(0.0, nrow(x), M))
    names(out) <- components_newdata$respVars
    for (i in seq_len(M)) {
        eta_i <-
            linpred_long(components_newdata$X, lapply(betas, i_row, i),
                         components_newdata$Z, splt_REs(b_mat[, , i], ind_RE),
                         components_newdata$id, level = level)
        for (j in seq_len(K)) {
            out[[j]][, i] <- if (pred_type == "response") {
                mu_fun(eta_i[[j]], links[j])
            } else eta_i[[j]]
        }
    }
    res1 <- list(preds = lapply(out, rowMeans, na.rm = TRUE),
                low = lapply(out, rowQuantile, probs = (1 - CI_level) / 2),
                upp = lapply(out, rowQuantile, probs = (1 + CI_level) / 2))
    if (return_newdata) {
        n <- nrow(newdata)
        preds <- mapply2(fix_NAs_preds, res1$preds, components_newdata$NAs,
                         MoreArgs = list(n = n))
        names(preds) <- paste0("pred_", components_newdata$respVars)
        low <- mapply2(fix_NAs_preds, res1$low, components_newdata$NAs,
                         MoreArgs = list(n = n))
        names(low) <- paste0("low_", components_newdata$respVars)
        upp <- mapply2(fix_NAs_preds, res1$upp, components_newdata$NAs,
                         MoreArgs = list(n = n))
        names(upp) <- paste0("upp_", components_newdata$respVars)
        l <- c(preds, low, upp)
        l <- l[c(matrix(seq_along(l), ncol = length(preds), byrow = TRUE))]
        res1 <- cbind(newdata, as.data.frame(do.call("cbind", l)))
    }
    ############################################################################
    ############################################################################
    # Predictions for newdata2
    if (!is.null(newdata2)) {
        terms_FE_noResp <- object$model_info$terms$terms_FE_noResp
        terms_RE <- object$model_info$terms$terms_RE
        mf_FE <- lapply(terms_FE_noResp, model.frame.default, data = newdata2)
        mf_RE <- lapply(terms_RE, model.frame.default, data = newdata2)
        NAs_FE <- lapply(mf_FE, attr, "na.action")
        NAs_RE <- lapply(mf_RE, attr, "na.action")
        mf_FE <- mapply2(fix_NAs_fixed, mf_FE, NAs_FE, NAs_RE)
        mf_RE <- mapply2(fix_NAs_random, mf_RE, NAs_RE, NAs_FE)
        X <- mapply2(model.matrix.default, terms_FE_noResp, mf_FE)
        Z <- mapply2(model.matrix.default, terms_RE, mf_RE)
        NAs <- mapply2(c, NAs_FE, NAs_RE)
        idL <- newdata2[[object$model_info$var_names$idVar]]
        unq_id <- unique(idL)
        idL <- mapply2(exclude_NAs, NAs_FE, NAs_RE, MoreArgs = list(id = idL))
        idL <- lapply(idL, match, table = unq_id)
        out <- lapply(X, function (x) matrix(0.0, nrow(x), M))
        names(out) <- components_newdata$respVars
        for (i in seq_len(M)) {
            eta_i <-
                linpred_long(X, lapply(betas, i_row, i), Z,
                             splt_REs(b_mat[, , i], ind_RE), idL, level = level)
            for (j in seq_len(K)) {
                out[[j]][, i] <- if (pred_type == "response") {
                    mu_fun(eta_i[[j]], links[j])
                } else eta_i[[j]]
            }
        }
        res2 <- list(preds = lapply(out, rowMeans, na.rm = TRUE),
                     low = lapply(out, rowQuantile, probs = (1 - CI_level) / 2),
                     upp = lapply(out, rowQuantile, probs = (1 + CI_level) / 2))
        if (return_newdata) {
            n <- nrow(newdata2)
            preds <- mapply2(fix_NAs_preds, res2$preds, NAs,
                             MoreArgs = list(n = n))
            names(preds) <- paste0("pred_", components_newdata$respVars)
            low <- mapply2(fix_NAs_preds, res2$low, NAs,
                           MoreArgs = list(n = n))
            names(low) <- paste0("low_", components_newdata$respVars)
            upp <- mapply2(fix_NAs_preds, res2$upp, NAs,
                           MoreArgs = list(n = n))
            names(upp) <- paste0("upp_", components_newdata$respVars)
            l <- c(preds, low, upp)
            l <- l[c(matrix(seq_along(l), ncol = length(preds), byrow = TRUE))]
            res2 <- cbind(newdata2, as.data.frame(do.call("cbind", l)))
        }
    }
}





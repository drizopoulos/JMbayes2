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
newdata <- pbc2[pbc2$id %in% c(2, 3, 10), ]
ind <- tapply(row.names(newdata), factor(newdata$id), tail, 1)
newdata2 <- newdata[ind[!is.na(ind)], ]; rm(ind)
newdata2 <- newdata2[rep(1:nrow(newdata2), each = 3), ]
rownames(newdata2) <- seq(1:nrow(newdata2))
newdata2$year <- with(newdata2, ave(year, id, FUN = function (x) x + seq_along(x)))
process <- c("Event")
type <- "fixed"
level <- 1L
CI_level <- 0.95
pred_type <- "response"

#############################################################
#############################################################

# control
control <- object$control

# check for tibbles
if (inherits(newdata, "tbl_df") || inherits(newdata, "tbl")) {
    newdata <- as.data.frame(newdata)
}

# extract idVar and time_var
idVar <- object$model_info$var_names$idVar
time_var <- object$model_info$var_names$time_var

# set dataL as newdata; almost the same code as in jm()
dataL <- newdata
idL <- dataL[[idVar]]
nY <- length(unique(idL))
# order data by idL and time_var
if (is.null(dataL[[time_var]])) {
    stop("the variable specified in agument 'time_var' cannot be found ",
         "in the database of the longitudinal models.")
}
dataL <- dataL[order(idL, dataL[[time_var]]), ]

# extract terms
respVars <- object$model_info$var_names$respVars
terms_FE <- object$model_info$terms$terms_FE
terms_FE_noResp <- object$model_info$terms$terms_FE_noResp
terms_RE <- object$model_info$terms$terms_RE
terms_Surv <- object$model_info$terms$terms_Surv_noResp
Xbar <- object$model_data$Xbar
# create model frames
mf_FE_dataL <- lapply(terms_FE, model.frame.default, data = dataL)
mf_RE_dataL <- lapply(terms_RE, model.frame.default, data = dataL)

# we need to account for missing data in the fixed and random effects model frames,
# in parallel across outcomes (i.e., we will allow that some subjects may have no data
# for some outcomes)
NAs_FE_dataL <- lapply(mf_FE_dataL, attr, "na.action")
NAs_RE_dataL <- lapply(mf_RE_dataL, attr, "na.action")
mf_FE_dataL <- mapply(fix_NAs_fixed, mf_FE_dataL, NAs_FE_dataL, NAs_RE_dataL,
                      SIMPLIFY = FALSE)
mf_RE_dataL <- mapply(fix_NAs_random, mf_RE_dataL, NAs_RE_dataL, NAs_FE_dataL,
                      SIMPLIFY = FALSE)

# create response vectors
y <- lapply(mf_FE_dataL, model.response)
y <- lapply(y, function (yy) {
    if (is.factor(yy)) as.numeric(yy != levels(yy)[1L]) else yy
})
y[] <- lapply(y, as.matrix)
families <- object$model_info$families
family_names <- sapply(families, "[[", "family")
links <- sapply(families, "[[", "link")
# for family = binomial and when y has two columns, set the second column
# to the number of trials instead the number of failures
binomial_data <- family_names %in% c("binomial", "beta binomial")
trials_fun <- function (y) {
    if (NCOL(y) == 2L) y[, 2L] <- y[, 1L] + y[, 2L]
    y
}
y[binomial_data] <- lapply(y[binomial_data], trials_fun)
unq_id <- unique(idL)
idL <- mapply(exclude_NAs, NAs_FE_dataL, NAs_RE_dataL,
              MoreArgs = list(id = idL), SIMPLIFY = FALSE)
idL <- lapply(idL, match, table = unq_id)
idL_lp <- lapply(idL, function (x) match(x, unique(x)))
unq_idL <- lapply(idL, unique)
X <- mapply(model.matrix.default, terms_FE, mf_FE_dataL, SIMPLIFY = FALSE)
Z <- mapply(model.matrix.default, terms_RE, mf_RE_dataL, SIMPLIFY = FALSE)

################################

# extract terms
terms_Surv <- object$model_info$terms$terms_Surv
terms_Surv_noResp <- object$model_info$terms$terms_Surv_noResp
type_censoring <- object$model_info$type_censoring
dataS <- newdata
idT <- dataS[[idVar]]
mf_surv_dataS <- model.frame.default(terms_Surv, data = dataS)
if (!is.null(NAs_surv <- attr(mf_surv_dataS, "na.action"))) {
    idT <- idT[-NAs_surv]
    dataS <- dataS[-NAs_surv, ]
}
idT <- factor(idT, levels = unique(idT))
nT <- length(unique(idT))
if (nY != nT) {
    stop("the number of groups/subjects in the longitudinal and survival datasets ",
         "do not seem to match. A potential reason why this may be happening is ",
         "missing data in some covariates used in the individual models.")
}
Surv_Response <- model.response(mf_surv_dataS)
if (type_censoring == "right") {
    Time_right <- unname(Surv_Response[, "time"])
    Time_left <- Time_start <- trunc_Time <- rep(0.0, nrow(dataS))
    delta <-  unname(Surv_Response[, "status"])
} else if (type_censoring == "counting") {
    Time_start <- unname(Surv_Response[, "start"])
    Time_stop <- unname(Surv_Response[, "stop"])
    delta <-  unname(Surv_Response[, "status"])
    Time_right <- Time_stop
    trunc_Time <- Time_start # possible left truncation time
    Time_left <- rep(0.0, nrow(dataS))
} else if (type_censoring == "interval") {
    Time1 <-  unname(Surv_Response[, "time1"])
    Time2 <-  unname(Surv_Response[, "time2"])
    trunc_Time <- Time_start <- rep(0.0, nrow(dataS))
    delta <- unname(Surv_Response[, "status"])
    Time_right <- Time1
    Time_right[delta == 3] <- Time2[delta == 3]
    Time_right[delta == 2] <- 0.0
    Time_left <- Time1
    Time_left[delta <= 1] <- 0.0
}
if (type_censoring != "counting") {
    names(Time_right) <- names(Time_left) <- names(Time_start) <- idT
}
which_event <- which(delta == 1)
which_right <- which(delta == 0)
which_left <- which(delta == 2)
which_interval <- which(delta == 3)
# extract strata if present otherwise all subjects in one stratum
ind_strata <- attr(terms_Surv, "specials")$strata
strata <- if (is.null(ind_strata)) {
    rep(1, nrow(mf_surv_dataS))
} else {
    unclass(mf_surv_dataS[[ind_strata]])
}
Time_integration <- Time_right
Time_integration[which_left] <- Time_left[which_left]
Time_integration[which_interval] <- Time_left[which_interval]
Time_integration2 <- rep(0.0, length(Time_integration))
if (length(which_interval)) {
    Time_integration2[which_interval] <- Time_right[which_interval]
}

# create Gauss Kronrod points and weights
GK <- gaussKronrod(control$GK_k)
sk <- GK$sk
P <- c(Time_integration - trunc_Time) / 2
st <- outer(P, sk) + (c(Time_integration + trunc_Time) / 2)
log_Pwk <- unname(rep(log(P), each = length(sk)) +
                      rep_len(log(GK$wk), length.out = length(st)))
if (length(which_interval)) {
    # we take the absolute value because for the subjects for whom we do not have
    # interval censoring P2 will be negative and this will produce a NA when we take
    # the log in 'log_Pwk2'
    P2 <- abs(Time_integration2 - Time_integration) / 2
    st2 <- outer(P2, sk) + (c(Time_integration2 + Time_integration) / 2)
    log_Pwk2 <- rep(log(P2), each = length(sk)) +
        rep_len(log(GK$wk), length.out = length(st2))
} else {
    P2 <- st2 <- log_Pwk2 <- rep(0.0, nT * control$GK_k)
}

# knots for the log baseline hazard function
knots <- control$knots

# indices
ni_event <- tapply(idT, idT, length)
ni_event <- cbind(c(0, head(cumsum(ni_event), -1)), cumsum(ni_event))
id_H <- rep(paste0(idT, "_", unlist(tapply(idT, idT, seq_along))),
            each = control$GK_k)
id_H <- match(id_H, unique(id_H))
# id_H_ repeats each unique idT the number of quadrature points
id_H_ <- rep(idT, each = control$GK_k)
id_H_ <- match(id_H_, unique(id_H_))
id_h <- unclass(idT)

# Functional forms
functional_forms <- object$model_info$functional_forms
FunForms_per_outcome <- object$model_info$FunForms_per_outcome
collapsed_functional_forms <- object$model_info$collapsed_functional_forms
FunForms_cpp <- object$model_info$FunForms_cpp
FunForms_ind <- object$model_info$FunForms_ind
Funs_FunForms <- object$model_info$Funs_FunForms
eps <- object$model_info$eps
direction <- object$model_info$direction

# Design matrices
strata_H <- rep(strata, each = control$GK_k)
W0_H <- create_W0(c(t(st)), knots, control$Bsplines_degree + 1, strata_H)
dataS_H <- SurvData_HazardModel(st, dataS, Time_start,
                                paste0(idT, "_", strata), time_var)
mf <- model.frame.default(terms_Surv_noResp, data = dataS_H)
W_H <- construct_Wmat(terms_Surv_noResp, mf)
any_gammas <- as.logical(ncol(W_H))
if (!any_gammas) {
    W_H <- matrix(0.0, nrow = nrow(W_H), ncol = 1L)
}
attr <- lapply(functional_forms, extract_attributes, data = dataS_H)
eps <- lapply(attr, "[[", 1L)
direction <- lapply(attr, "[[", 2L)
X_H <- design_matrices_functional_forms(st, terms_FE_noResp,
                                        dataL, time_var, idVar, idT,
                                        collapsed_functional_forms, Xbar,
                                        eps, direction)
Z_H <- design_matrices_functional_forms(st, terms_RE,
                                        dataL, time_var, idVar, idT,
                                        collapsed_functional_forms, NULL,
                                        eps, direction)
U_H <- lapply(functional_forms, construct_Umat, dataS = dataS_H)
if (length(which_event)) {
    W0_h <- create_W0(Time_right, knots, control$Bsplines_degree + 1,
                      strata)
    dataS_h <- SurvData_HazardModel(Time_right, dataS, Time_start,
                                    paste0(idT, "_", strata), time_var)
    mf <- model.frame.default(terms_Surv_noResp, data = dataS_h)
    W_h <- construct_Wmat(terms_Surv_noResp, mf)
    if (!any_gammas) {
        W_h <- matrix(0.0, nrow = nrow(W_h), ncol = 1L)
    }
    X_h <- design_matrices_functional_forms(Time_right, terms_FE_noResp,
                                            dataL, time_var, idVar, idT,
                                            collapsed_functional_forms, Xbar,
                                            eps, direction)
    Z_h <- design_matrices_functional_forms(Time_right, terms_RE,
                                            dataL, time_var, idVar, idT,
                                            collapsed_functional_forms, NULL,
                                            eps, direction)
    U_h <- lapply(functional_forms, construct_Umat, dataS = dataS_h)
} else {
    W0_h <- W_h <- matrix(0.0)
    X_h <- Z_h <- U_h <- rep(list(matrix(0.0)), length(respVars))
}
if (length(which_interval)) {
    W0_H2 <- create_W0(c(t(st2)), knots, control$Bsplines_degree + 1,
                       strata_H)
    dataS_H2 <- SurvData_HazardModel(st2, dataS, Time_start,
                                     paste0(idT, "_", strata), time_var)
    mf2 <- model.frame.default(terms_Surv_noResp, data = dataS_H2)
    W_h <- construct_Wmat(terms_Surv_noResp, mf2)
    if (!any_gammas) {
        W_H2 <- matrix(0.0, nrow = nrow(W_H2), ncol = 1L)
    }
    X_H2 <- design_matrices_functional_forms(st, terms_FE_noResp,
                                             dataL, time_var, idVar, idT,
                                             collapsed_functional_forms, Xbar,
                                             eps, direction)
    Z_H2 <- design_matrices_functional_forms(st, terms_RE,
                                             dataL, time_var, idVar, idT,
                                             collapsed_functional_forms, NULL,
                                             eps, direction)
    U_H <- lapply(functional_forms, construct_Umat, dataS = dataS_H2)
} else {
    W0_H2 <- W_H2 <- matrix(0.0)
    X_H2 <- Z_H2 <- U_H2 <- rep(list(matrix(0.0)), length(respVars))
}
X_H[] <- lapply(X_H, docall_cbind)
X_h[] <- lapply(X_h, docall_cbind)
X_H2[] <- lapply(X_H2, docall_cbind)
Z_H[] <- lapply(Z_H, docall_cbind)
Z_h[] <- lapply(Z_h, docall_cbind)
Z_H2[] <- lapply(Z_H2, docall_cbind)


# MCMC sample
b <- lapply(sapply(Z, ncol), function (nc) matrix(0.0, nY, nc))
M <- sum(sapply(object$mcmc$bs_gammas, nrow))
get_param <- function (nam) {
    tht <- object$mcmc[[nam]]
    if (!is.null(tht)) docall_rbind(tht) else matrix(0.0, M, 1)
}
ind_betas <- grep("^betas", names(object$mcmc))
mcmc <- list(
    b = b, bs_gammas = get_param("bs_gammas"), gammas = get_param("gammas"),
    alphas = get_param("alphas"),
    betas = lapply(object$mcmc[ind_betas], docall_rbind)
)
has_sigmas <- object$model_data$has_sigmas
mcmc$sigmas <- matrix(0.0, M, length(has_sigmas))
mcmc$sigmas[, has_sigmas > 0] <- get_param("sigmas")
D <- get_param("D")
mcmc$D <- array(0.0, c(dim(lowertri2mat(D[1L, ])), M))
for (i in seq_len(M)) {
    mcmc$D[, , i] <- lowertri2mat(D[i, ])
}

# Priors
priors <- object$priors


Data <- list(
    ind_RE = object$model_data$ind_RE,
    X_H = X_H, X_h = X_h, X_H2 = X_H2,
    Z_H = Z_H, Z_h = Z_h, Z_H2 = Z_H2,
    U_H = U_H, U_h = U_h, U_H2 = U_H2,
    Wlong_bar = object$Wlong_bar, Wlong_sds = object$Wlong_sds,
    idT = match(idT, unique(idT)), log_Pwk = log_Pwk, log_Pwk2 = log_Pwk2,
    id_H = id_H, id_H_ = id_H_, id_h = id_h,
    which_event = which_event, which_right = which_right,
    which_left = which_left, which_interval = which_interval,
    ni_event = ni_event, FunForms_cpp = FunForms_cpp,
    FunForms_ind = FunForms_ind, Funs_FunForms = Funs_FunForms,
    X = X, Z = Z, y = y, family_names = family_names,
    links = links, extra_parms = object$model_data$extra_parms,
    unq_idL = unq_idL, idL_lp = idL_lp, idL = idL
)

control <- list(GK_k = object$control$GK_k)








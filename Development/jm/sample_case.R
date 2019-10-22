##########################################################################################
# Author: D. Rizopoulos                                                                  #
# Aim: Create a generic setting based on the PBC dataset to serve as an example          #
##########################################################################################

# load packages and data
library("survival")
library("nlme")
library("GLMMadaptive")
library("splines")
library("Formula")
data("pbc2", package = "JM")
data("pbc2.id", package = "JM")
source(file.path(getwd(), "Development/jm/help_functions.R"))

####################

# create artificial interval censored data
pbc2$status2 <- as.numeric(pbc2$status != "alive")
pbc2.id$status2 <- as.numeric(pbc2.id$status != "alive")
pbc2$status3 <- as.character(pbc2$status)
ff <- function (x) {
    out <- if (x[1L] %in% c('dead', 'transplanted')) 'dead' else
        switch(sample(1:3, 1), '1' = "alive", '2' = "left", '3' = "interval")
    rep(out, length(x))
}
pbc2$status3 <- unlist(with(pbc2, lapply(split(status3, id), ff)), use.names = FALSE)
pbc2$status3 <- unname(with(pbc2, sapply(status3, function (x)
    switch(x, 'dead' = 1, 'alive' = 0, 'left' = 2, 'interval' = 3))))
pbc2$yearsU <- as.numeric(NA)
pbc2$yearsU[pbc2$status3 == 3] <- pbc2$years[pbc2$status3 == 3] +
    runif(sum(pbc2$status3 == 3), 0, 4)
pbc2.id <- pbc2[!duplicated(pbc2$id), ]


fm1 <- lme(log(serBilir) ~ ns(year, 3, B = c(0, 11)) * sex + age, data = pbc2,
           random = ~ ns(year, 2, B = c(0, 11)) | id)
fm2 <- lme(serChol ~ year + sex + age, data = pbc2, random = ~ year | id,
           na.action = na.exclude)
fm3 <- mixed_model(hepatomegaly ~ year + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())
fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())

CoxFit <- coxph(Surv(years, status2) ~ ns(age, 3) + sex + cluster(id),
                data = pbc2.id, model = TRUE)

survFit <- survreg(Surv(years, yearsU, status3, type = "interval") ~ drug + age + cluster(id),
                   data = pbc2.id, model = TRUE)

##########################################################################################

# the arguments of the jm() function

Surv_object = survFit#CoxFit
Mixed_objects = list(fm1, fm2, fm3, fm4)
data_Surv = NULL
timeVar = "year"

# default function_form
functional_form = Formula(~ value(log(serBilir)) + slope(log(serBilir)) +
                              value(log(serBilir)):sex | value(serChol) | value(hepatomegaly)
                          | value(ascites))

# complex example function_form
functional_form = Formula(~ value(log(serBilir)) + slope(log(serBilir)) |
                              value(serChol) + value(serChol):sex |
                              logit(value(hepatomegaly)) |
                              value(ascites) + area(ascites))

##########################################################################################

con <- list(GK_k = 15L, Bsplines_degree = 2, base_hazard_segments = 10)


# extract the data from each of the mixed models
# and check whether the same data have been used;
# otherwise an error
datas <- lapply(Mixed_objects, "[[", "data")
if (!all(sapply(datas[-1L], all.equal, datas[[1L]]))) {
    stop("It seems that some of the mixed models have been fitted to different versions ",
         "of the dataset. Use the same exact dataset in the calls to lme() ",
         " and mixed_model().")
}
dataL <- datas[[1L]]
rm(datas)

# extract id variable (again we assume a single grouping variable)
id_names <- sapply(Mixed_objects, function (object)
    names(if (inherits(object, "MixMod")) object$id[1L] else object$groups[1L]))
if (!all(id_names == id_names[1L])) {
    stop("it seems that different grouping variables have been used in the mixed models.")
}
idVar <- id_names[1L]
idL <- dataL[[idVar]]
nY <- length(unique(idL))

# order data by idL
dataL <- dataL[order(idL), ]

# extract terms from mixed models
# (function extract_terms() is defined in help_functions)
terms_FE <- lapply(Mixed_objects, extract_terms, which = "fixed", data = dataL)
terms_FE_noResp <- lapply(terms_FE, delete.response)
terms_RE <- lapply(Mixed_objects, extract_terms, which = "random", data = dataL)

# create model frames
mf_FE_dataL <- lapply(terms_FE, model.frame.default, data = dataL)
mf_RE_dataL <- lapply(terms_RE, model.frame.default, data = dataL)

# we need to account for missing data in the fixed and random effects model frames,
# in parallel across outcomes (i.e., we will allow that some subjects may have no data
# for some outcomes)
NAs_FE_dataL <- lapply(mf_FE_dataL, attr, "na.action")
NAs_RE_dataL <- lapply(mf_RE_dataL, attr, "na.action")
mf_FE_dataL <- mapply(fix_NAs_fixed, mf_FE_dataL, NAs_FE_dataL, NAs_RE_dataL)
mf_RE_dataL <- mapply(fix_NAs_random, mf_RE_dataL, NAs_RE_dataL, NAs_FE_dataL)

# create response vectors
y <- lapply(mf_FE_dataL, model.response)

# exctract families
families <- lapply(Mixed_objects, "[[", "family")
families[sapply(families, is.null)] <- rep(list(gaussian()),
                                           sum(sapply(families, is.null)))
# create the idL per outcome
# IMPORTANT: some ids may be missing when some subjects have no data for a particular outcome
# This needs to be taken into account when using idL for indexing. Namely, a new id variable
# will need to be created in jm_fit()
unq_id <- unique(idL)
idL <- mapply(exclude_NAs, NAs_FE_dataL, NAs_RE_dataL,
             MoreArgs = list(id = idL), SIMPLIFY = FALSE)
idL <- lapply(idL, match, table = unq_id)
idL_lp <- lapply(idL, function (x) match(x, unique(x)))
unq_idL <- lapply(idL, unique)
sample_size_Long <- sapply(idL, function (x) length(unique(x)))

# create design matrices for mixed models
X <- mapply(model.matrix.default, terms_FE, mf_FE_dataL)
Z <- mapply(model.matrix.default, terms_RE, mf_RE_dataL)

########################################################

# We require users to include the id variable as cluster, even in the case of simple
# right censoring. The estimated coefficients are in either case the same. This will
# give us the id variable to match with the longitudinal data
if (is.null(Surv_object$model$cluster)) {
    stop("you need to refit the Cox or AFT model and include in the right hand side of the ",
         "formula the 'cluster()' function using as its argument the subjects' ",
         "id indicator. These ids need to be the same as the ones used to fit ",
         "the mixed effects model.\n")
}

# try to recover survival dataset
if (is.null(data_Surv))
    try(dataS <- eval(Surv_object$call$data, envir = parent.frame()),
        silent = TRUE)
if (inherits(dataS, "try-error")) {
    stop("could not recover the dataset used to fit the Cox model; please provide this ",
         "dataset in the 'data_Surv' argument of jm().")
}

# terms for survival model
terms_Surv <- Surv_object$terms
terms_Surv <- drop.terms(terms_Surv, attr(terms_Surv,"specials")$cluster - 1,
                         keep.response = TRUE)
terms_Surv_noResp <- delete.response(terms_Surv)
mf_surv_dataS <- model.frame.default(terms_Surv, data = dataS)

# survival times
Surv_Response <- model.response(mf_surv_dataS)
type_censoring <- attr(Surv_Response, "type")
idT <- Surv_object$model$cluster
idT <- factor(idT, levels = unique(idT))
if (!all(idT %in% dataL[[idVar]])) {
    stop("it seems that some of the levels of the 'cluster()' variable in the survival ",
         "object can be found in the dataset of the Mixed_objects. Please check that ",
         "the same subjects / groups are used in the datasets used to fit the mixed ",
         "and survival models.")
}
# we need to check that the ordering of the subjects in the same in dataL and dataS.
# If not, then a warning and do it internally
if (!all(order(unique(idT)) == order(unique(dataL[[idVar]])))) {
    warning("It seems that the ordering of the subjects in dataset used to fit the ",
            "mixed models and the dataset used for the survival model is not the same. ",
            "We set internally the datasets in the same order, but it would be best ",
            "that you do it beforehand on your own.")
    dataS <- dataS[order(idT), ]
    mf_surv_dataS <- model.frame.default(terms_Surv, data = dataS)
    Surv_Response <- model.response(mf_surv_dataS)
}
nT <- length(unique(idT))
if (nY != nT) {
    stop("the number of groups/subjects in the longitudinal and survival datasets ",
         "do not seem to match.")
}
# Notation:
#  - Time_right: event or right censoring time
#  - Time_left: left censoring time
#  - trunc_Time: truncation time
#  - delta: 0 of right censored, 1 for event, 2 for left censored,
#           3 for interval censored
if (type_censoring == "right") {
    Time_right <- unname(Surv_Response[, "time"])
    Time_left <- rep(0.0, nT)
    trunc_Time <- rep(0.0, nT)
    delta <-  unname(Surv_Response[, "status"])
} else if (type_censoring == "counting") {
    Time_start <-  unname(Surv_Response[, "start"])
    Time_stop <-  unname(Surv_Response[, "stop"])
    delta <-  unname(Surv_Response[, "status"])
    Time_right <- tapply(Time_stop, idT, tail, n = 1) # time of event
    trunc_Time <- tapply(Time_start, idT, head, n = 1) # possible left truncation time
    Time_left <- rep(0.0, nT)
    delta <- tapply(event, idT, tail, n = 1) # event indicator at Time_right
} else if (type_censoring == "interval") {
    Time1 <-  unname(Surv_Response[, "time1"])
    Time2 <-  unname(Surv_Response[, "time2"])
    trunc_Time <- rep(0.0, nT)
    delta <-  unname(Surv_Response[, "status"])
    Time_right <- Time1
    Time_right[delta == 3] <- Time2[delta == 3]
    Time_right[delta == 2] <- 0.0
    Time_left <- Time1
    Time_left[delta <= 1] <- 0.0
}
which_event <- which(delta == 1)
which_right <- which(delta == 0)
which_left <- which(delta == 2)
which_interval <- which(delta == 3)

# 'Time_integration' is the upper limit of the integral in likelihood
# of the survival model. For subjects with event (delta = 1), for subjects with
# right censoring and for subjects with interval censoring we need to integrate
# up to 'Time_right'. For subjects with left censoring we need to integrate up to
# 'Time_left'; hence we set for them 'Time_integration = Time_left'. For subjects
# with interval censoring we need also 'Time_integration2' which is equal to
# 'Time_left'
Time_integration <- Time_right
Time_integration[which_left] <- Time_left[which_left]
Time_integration2 <- rep(0.0, nT)
if (length(which_interval)) {
    Time_integration2[which_interval] <- Time_left[which_interval]
}

# create Gauss Kronrod points and weights
GK <- gaussKronrod(15)
wk <- GK$wk
sk <- GK$sk
P <- c(Time_integration - trunc_Time) / 2
st <- outer(P, sk) + (c(Time_integration + trunc_Time) / 2)
if (length(which_interval)) {
    P2 <- c(Time_integration2 - trunc_Time) / 2
    st2 <- outer(P2, sk) + (c(Time_integration2 + trunc_Time) / 2)
} else {
    P2 <- st2 <- NULL
}

# knots for the log baseline hazard function
if (is.null(con$knots)) {
    con$knots <- knots(0, floor(max(Time_integration)) + 1,
                       con$base_hazard_segments, con$Bsplines_degree)
}

# Extract functional forms per longitudinal outcome
if (length(functional_form)[2L] != length(Mixed_objects)) {
    stop("it seems that you have not specified any functional form for some of the ",
         "longitudinal outcomes.")
}
long_resp_vars <- sapply(Mixed_objects,
                         function (obj) as.character(formula(terms(obj)))[2L])
long_resp_var_in_functional_form <- sapply(seq_along(Mixed_objects),
        function (i) as.character(formula(functional_form, rhs = i))[2L])
ordering_of_outcomes <- sapply(long_resp_vars, grep, x = long_resp_var_in_functional_form,
                               fixed = TRUE)
functional_forms_per_outcome <- lapply(ordering_of_outcomes,
                                       extract_functional_forms_per_outcome)
collapsed_functional_forms <- lapply(functional_forms_per_outcome,
                                     function (x) names(x[sapply(x, length) > 0]))

# design matrices for the survival submodel:
#  - W0 is the design matrix for the log baseline hazard
#  - W is the design matrix for the covariates in the Surv_object
#    (including time-varying covariates)
#  - X is the design matrix for the fixed effects, per outcome and functional form
#  - Z is the design matrix for the random effects, per outcome and functional form
# in the above design matrices we put the "_h" to denote calculation at the event time
# 'Time_right', we put "_H" to denote calculation at the 'Time_integration', and
# "_H2" to denote calculation at the 'Time_integration2'.
if (length(which_event)) {
    W0_h <- splineDesign(con$knots, Time_right, ord = con$Bsplines_degree + 1)
    W_h <- model.matrix.default(terms_Surv_noResp, mf_surv_dataS)[, -1, drop = FALSE]
    X_h <- desing_matrices_functional_forms(Time_right, terms_FE_noResp,
                                            dataL, timeVar, idVar,
                                            collapsed_functional_forms)
    Z_h <- desing_matrices_functional_forms(Time_right, terms_RE,
                                            dataL, timeVar, idVar,
                                            collapsed_functional_forms)
}
W0_H <- splineDesign(con$knots, c(t(st)), ord = con$Bsplines_degree + 1)
dataS_int <- SurvData_HazardModel(st, dataS, Time_start, idT) # <- check that it works with start stop Surv
mf <- model.frame.default(terms_Surv_noResp, data = dataS_int)
W_H <- model.matrix.default(terms_Surv_noResp, mf)[, -1, drop = FALSE]
X_H <- desing_matrices_functional_forms(st, terms_FE_noResp,
                                        dataL, timeVar, idVar,
                                        collapsed_functional_forms)
Z_H <- desing_matrices_functional_forms(st, terms_RE,
                                        dataL, timeVar, idVar,
                                        collapsed_functional_forms)
if (length(which_interval)) {
    W0_H2 <- splineDesign(con$knots, c(t(st2)), ord = con$Bsplines_degree + 1)
    dataS_int2 <- SurvData_HazardModel(st2, dataS, Time_start, idT)
    mf2 <- model.frame.default(terms_Surv_noResp, data = dataS_int2)
    W_H2 <- model.matrix.default(terms_Surv_noResp, mf2)[, -1, drop = FALSE]
    X_H2 <- desing_matrices_functional_forms(st, terms_FE_noResp,
                                             dataL, timeVar, idVar,
                                             collapsed_functional_forms)
    Z_H2 <- desing_matrices_functional_forms(st, terms_RE,
                                             dataL, timeVar, idVar,
                                             collapsed_functional_forms)
}

#############################################################

# extract initial values
betas <- lapply(Mixed_objects, fixef)
D <- bdiag(lapply(Mixed_objects, extract_D))
b <- mapply(extract_b, Mixed_objects, unq_idL, MoreArgs = list(n = nY),
            SIMPLIFY = FALSE)
gammas <- coef(Surv_object)

################################################################################

linpred_mixed <- function (X, betas, Z, b, id) {
    n_outcomes <- length(X)
    out <- vector("list", n_outcomes)
    for (i in seq_len(n_outcomes)) {
        X_i <- X[[i]]
        betas_i <- betas[[i]]
        Z_i <- Z[[i]]
        b_i <- b[[i]]
        id_i <- id[[i]]
        out[[i]] <- X_i %*% betas_i + rowSums(Z_i * b_i[id_i, ])
    }
    out
}

xxx <- linpred_mixed(X, betas, Z, b, idL_lp)

log_dens_Funs <- function (family) {
    gaussian_log_dens <- function (y, eta, mu_fun = NULL, phis, eta_zi = NULL) {
        dnorm(y, eta, exp(phis), log = TRUE)
    }
    switch(family$family,
           "gaussian" = gaussian_log_dens,
           "binomial" = GLMMadaptive:::binomial_log_dens,
           "poisson" = GLMMadaptive:::poisson_log_dens)
}

lapply(families, log_dens_Funs)

log_density_mixed <- function (y, linear_predictor, sigmas, Funs, n, id) {
    n_outcomes <- length(y)
    out <- numeric(n)
    for (i in seq_len(n_outcomes)) {
        y_i <- y[[i]]
        eta_i <- linear_predictor[[i]]
        sigma_i <- sigmas[[i]]
        id_i <- id[[i]]
        log_dens_i <- Funs[[i]]
        out[id_i] <- out[id_i] + rowsum(log_dens_i(y_i, eta_i, sigma_i), id_i) # <- fix the second id
    }
    out
}






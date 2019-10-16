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

Surv_object = CoxFit
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

con <- list("slope_eps" = 0.001, GK_k = 15L)


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
# Notation:
#  - Time_right: event or right censoring time
#  - Time_left: left censoring time
#  - trunc_Time: truncation time
#  - delta: 0 of right censored, 1 for event, 2 for left censored,
#           3 for interval censored
if (type_censoring == "right") {
    Time_right <- Surv_Response[, "time"]
    Time_left <- rep(as.numeric(NA), nT)
    trunc_Time <- rep(0.0, nT)
    delta <- Surv_Response[, "status"]
} else if (type_censoring == "counting") {
    Time_start <- Surv_Response[, "start"]
    Time_stop <- Surv_Response[, "stop"]
    event <- Surv_Response[, "status"]
    Time_right <- tapply(Time_stop, idT, tail, n = 1) # time of event
    trunc_Time <- tapply(Time_start, idT, head, n = 1) # possible left truncation time
    Time_left <- rep(as.numeric(NA), nT)
    event <- tapply(event, idT, tail, n = 1) # event indicator at Time_right
} else if (type_censoring == "interval") {
    Time1 <- Surv_Response[, "time1"]
    Time2 <- Surv_Response[, "time2"]
    trunc_Time <- rep(0.0, nT)
    event <- Surv_Response[, "status"]
    Time_right <- Time1
    Time_right[event == 3] <- Time2[event == 3]
    Time_right[event == 2] <- as.numeric(NA)
    Time_left <- Time1
    Time_left[event %in% c(0, 1)] <- as.numeric(NA)
}

# covariates design matrix Cox model
W <- model.matrix.default(terms_Surv, mf_surv_dataS)[, -1, drop = FALSE]


#############################################################

if (length(functional_form)[2L] != length(Mixed_objects)) {
    stop("it seems that you have not specified any functional form for some of the ",
         "longitudinal outcomes.")
}

long_resp_vars <- sapply(Mixed_objects,
                         function (obj) as.character(formula(terms(obj)))[2L])

Form <- formula(functional_form, rhs = 1)
long_resp_var_in_functional_form <- sapply(seq_along(Mixed_objects),
       function (i) as.character(formula(functional_form, rhs = i))[2L])
ordering_of_outcomes <- sapply(long_resp_vars, grep, x = long_resp_var_in_functional_form,
                               fixed = TRUE)

functional_forms_per_outcome <- lapply(ordering_of_outcomes,
                                       extract_functional_forms_per_outcome)
collapsed_functional_forms <- lapply(functional_forms_per_outcome,
                                     function (x) names(x[sapply(x, length) > 0]))

# extract the last row per subject, and set the timeVar equal to Time

desing_matrices_functional_forms <- function (time, dataL, timeVar, idVar,
                                              terms_FE_noResp, terms_RE) {
    desgn_matr <- function (time, terms) {
        D <- LongData_HazardModel(time, dataL, dataL[[timeVar]],
                                  dataL[[idVar]], timeVar)
        mf <- lapply(terms, model.frame.default, data = D)
        mapply(model.matrix.default, terms, mf)
    }
    degn_matr_slp <- function (time, terms) {
        M1 <- desgn_matr(time + 0.001, terms)
        M2 <- desgn_matr(time - 0.001, terms)
        mapply(function (x1, x2) (x1 - x2) / 0.002, M1, M2)
    }
    degn_matr_area <- function (time, terms) {
        GK <- gaussKronrod(15)
        wk <- GK$wk
        sk <- GK$sk
        P <- time / 2
        st <- outer(P, sk + 1)
        out <- vector("list", 15L)
        for (i in seq_len(15)) {
            M <- desgn_matr(st[, i], terms)
            out[[i]] <- lapply(M, "*", wk[i])
        }
        out <- P * Reduce("+", out)  # <------------------ FIX not correct
    }
    ################
    X_value <- desgn_matr(time, terms_FE_noResp)
    Z_value <- desgn_matr(time, terms_RE)
    ################
    X_slope <- degn_matr_slp(time, terms_FE_noResp)
    Z_slope <- degn_matr_slp(time, terms_RE)
    ################
    X_area <- degn_matr_area(time, terms_FE_noResp)
    Z_area <- degn_matr_area(time, terms_RE)
}

dataL_Time_right <- LongData_HazardModel(Time_right, dataL, dataL[[timeVar]],
                                         dataL[[idVar]], timeVar)
mf_FE_dataXX <- lapply(terms_FE_noResp, model.frame.default, data = dataXX)
X_h <- mapply(model.matrix.default, terms_FE_noResp, mf_FE_dataXX)




dataL_at_Time_plus_eps <- LongData_HazardModel(Time + con$slope_eps,
                                               dataL, dataL[[timeVar]], dataL[[idVar]],
                                               timeVar)
dataL_at_Time_min_eps <- LongData_HazardModel(Time - con$slope_eps,
                                               dataL, dataL[[timeVar]], dataL[[idVar]],
                                               timeVar)

xxx <- LongData_HazardModel(Time, dataL, dataL[[timeVar]], dataL[[idVar]],
                            timeVar)

all.equal(xxx, dataL_at_Time)

#############################################################


# check if the max(sample_size_Long) == sample size surv

# extract initial values
betas <- lapply(Mixed_objects, fixef)



################################################################################

GQsurv <- gaussKronrod(con$GK_k)
wk <- GQsurv$wk
sk <- GQsurv$sk
K <- length(sk)
P <- (Time_right - trunc_Time) / 2
GK_times <- outer(P, sk) + c(Time_right + trunc_Time) / 2
GK_weights <- rep(P, each = con$GK_k) * rep(wk, nT)

dataXX <- LongData_HazardModel(GK_times, dataL, dataL[[timeVar]],
                                dataL[[idVar]], timeVar)

terms_FE_noResp <- lapply(terms_FE, delete.response)
mf_FE_dataXX <- lapply(terms_FE_noResp, model.frame.default, data = dataXX)
XX <- mapply(model.matrix.default, terms_FE_noResp, mf_FE_dataXX)

xxx <- c(XX[[1]] %*% betas[[1]])
sum(GK_weights * exp(xxx))

idGK <- rep(1:nT, each = 15)
sum(exp(rowsum(log(GK_weights) + XX[[1]], idGK) %*% betas[[1]]))



##########################################################################################
# Author: D. Rizopoulos                                                                  #
# Aim: Create a generic setting based on the PBC dataset to serve as an example          #
##########################################################################################

# load packages and data
library("survival")
library("nlme")
library("GLMMadaptive")
library("splines")
data("pbc2", package = "JM")
data("pbc2.id", package = "JM")
source(file.path(getwd(), "Development/jm/help_functions.R"))


fm1 <- lme(log(serBilir) ~ ns(year, 3, B = c(0, 11)) * sex + age, data = pbc2,
           random = ~ ns(year, 2, B = c(0, 11)) | id)
fm2 <- lme(log(serChol) ~ year + sex + age, data = pbc2, random = ~ year | id,
           na.action = na.exclude)
fm3 <- mixed_model(hepatomegaly ~ year + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())
fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())

CoxFit <- coxph(Surv(years, status2) ~ ns(age, 3) + sex + cluster(id),
                data = pbc2.id, model = TRUE)

##########################################################################################

# the arguments of the jm() function

Cox_object = CoxFit
Mixed_objects = list(fm1, fm2, fm3, fm4)
data_Surv = NULL
timeVar = "year"

##########################################################################################

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
if (is.null(Cox_object$model$cluster)) {
    stop("you need to refit the Cox and include in the right hand side of the ",
         "formula the 'cluster()' function using as its argument the subjects' ",
         "id indicator. These ids need to be the same as the ones used to fit ",
         "the mixed effects model.\n")
}

# try to recover survival dataset
if (is.null(data_Surv))
    try(dataS <- eval(Cox_object$call$data, envir = parent.frame()),
        silent = TRUE)
if (inherits(dataS, "try-error")) {
    stop("could not recover the dataset used to fit the Cox model; please provide this ",
         "dataset in the 'data_Surv' argument of jm().")
}

# terms for survival model
terms_Surv <- Cox_object$terms
mf_surv_dataS <- model.frame.default(terms_Surv, data = dataS)

# survival times
Surv_Response <- model.response(mf_surv_dataS)
type_censoring <- attr(Surv_Response, "type")


# covariates design matrix Cox model
W <- model.matrix.default(terms_Surv, mf_surv_dataS)[, -1, drop = FALSE]

# check if the max(sample_size_Long) == sample size surv

# extract initial values
betas <- lapply(Mixed_objects, fixef)











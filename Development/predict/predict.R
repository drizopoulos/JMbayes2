if (FALSE) {
    library("JMbayes2")
    pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
    CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
    fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
    fm2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)
    fm3 <- mixed_model(ascites ~ year + sex, data = pbc2,
                       random = ~ year | id, family = binomial())

    jointFit1 <- jm(CoxFit, list(fm1, fm2, fm3), time_var = "year")

    fix_NAs_fixed <- JMbayes2:::fix_NAs_fixed
    fix_NAs_random <- JMbayes2:::fix_NAs_random
    exclude_NAs <- JMbayes2:::exclude_NAs
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
terms_FE <- object$model_info$terms$terms_FE
terms_RE <- object$model_info$terms$terms_RE
terms_Surv <- object$model_info$terms$terms_Surv_noResp
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
if (!all(idT %in% dataL[[idVar]])) {
    stop("it seems that some of the levels of the id variable in the survival dataset",
         "cannot be found in the dataset of the Mixed_objects. Please check that ",
         "the same subjects/groups are used in the datasets used to fit the mixed ",
         "and survival models. Also, the name of the subjects/groups variable ",
         "in the different datasets used to fit the individual models ",
         "needs to be the same in all of the datasets.")
}
# we need to check that the ordering of the subjects in the same in dataL and dataS.
# If not, then a warning and do it internally
if (!all(order(unique(idT)) == order(unique(dataL[[idVar]])))) {
    warning("It seems that the ordering of the subjects in the dataset used to fit the ",
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
         "do not seem to match. A potential reason why this may be happening is ",
         "missing data in some covariates used in the individual models.")
}






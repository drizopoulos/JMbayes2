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

# extract terms
terms_FE <- object$model_info$terms$terms_FE
terms_RE <- object$model_info$terms$terms_RE
terms_Surv <- object$model_info$terms$terms_Surv_noResp
# create model frames
mf_FE_dataL <- lapply(terms_FE, model.frame.default, data = newdata)
mf_RE_dataL <- lapply(terms_RE, model.frame.default, data = newdata)

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
binomial_data <- family_names == "binomial"
trials_fun <- function (y) {
    if (NCOL(y) == 2L) y[, 2L] <- y[, 1L] + y[, 2L]
    y
}
y[binomial_data] <- lapply(y[binomial_data], trials_fun)
idVar <- object$model_info$var_names$idVar
idL <- newdata[[idVar]] # error if not present
unq_id <- unique(idL)
idL <- mapply(exclude_NAs, NAs_FE_dataL, NAs_RE_dataL,
              MoreArgs = list(id = idL), SIMPLIFY = FALSE)
idL <- lapply(idL, match, table = unq_id)
idL_lp <- lapply(idL, function (x) match(x, unique(x)))
unq_idL <- lapply(idL, unique)
X <- mapply(model.matrix.default, terms_FE, mf_FE_dataL, SIMPLIFY = FALSE)
Z <- mapply(model.matrix.default, terms_RE, mf_RE_dataL, SIMPLIFY = FALSE)







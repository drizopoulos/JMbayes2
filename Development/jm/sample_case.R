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

CoxFit <- coxph(Surv(years, status2) ~ ns(age, 3) + sex, data = pbc2.id, model = TRUE)

##########################################################################################

# the arguments of the jm() function


Cox_object = CoxFit
Mixed_objects = list(fm1, fm2, fm3, fm4)


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
data <- datas[[1L]]
rm(datas)

# extract terms from mixed models
# (function extract_terms() is defined in help_functions)
terms_FE <- lapply(Mixed_objects, extract_terms, which = "fixed")
terms_RE <- lapply(Mixed_objects, extract_terms, which = "random")

# create model frames
mf_FE_data <- lapply(terms_FE, model.frame.default, data = data)
mf_RE_data <- lapply(terms_RE, model.frame.default, data = data)

# we need to account for missing data in the fixed and random effects model frames,
# in parallel across outcomes (i.e., we will allow that some subjects may have no data
# for some outcomes)
NAs_FE_data <- lapply(mf_FE_data, attr, "na.action")
NAs_RE_data <- lapply(mf_RE_data, attr, "na.action")
mf_FE_data <- mapply(fix_NAs_fixed, mf_FE_data, NAs_FE_data, NAs_RE_data)
mf_RE_data <- mapply(fix_NAs_random, mf_RE_data, NAs_RE_data, NAs_FE_data)

# create response vectors
y <- lapply(mf_FE_data, model.response)

# exctract families
families <- lapply(Mixed_objects, "[[", "family")
families[sapply(families, is.null)] <- rep(list(gaussian()), sum(sapply(families, is.null)))

# create design matrices for mixed models
X <- mapply(model.matrix.default, terms_FE, mf_FE_data)
Z <- mapply(model.matrix.default, terms_RE, mf_RE_data)














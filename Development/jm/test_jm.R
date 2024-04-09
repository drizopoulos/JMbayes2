##########################################################################################
# Author: D. Rizopoulos                                                                  #
# Aim: Test function jm()                                                                #
#########################################################################################
library("survival")
library("nlme")
library("GLMMadaptive")
library("coda")
library("lattice")
library("splines")
library("matrixStats")
source("./R/jm.R")
source("./R/jm_fit.R")
source("./R/help_functions.R")
source("./R/basic_methods.R")
source("./R/predict_funs.R")
source("./R/create_Wlong_mats.R")
source("./Development/jm/PBC_data.R")
Rcpp::sourceCpp('src/mcmc_fit.cpp')


#
pbc2$prothrombin[pbc2$id == levels(pbc2$id)[1L]] <- NA
pbc2$prothrombin[pbc2$id == levels(pbc2$id)[2L]] <- NA

fm1 <- lme(log(serBilir) ~ ns(year, 2, B = c(0, 15)) * sex, data = pbc2,
           random = list(id = pdDiag(form = ~ ns(year, 2, B = c(0, 15)))))

fm2 <- lme(prothrombin ~ ns(year, 2, B = c(0, 15)) * sex, data = pbc2,
           random = ~ ns(year, 2, B = c(0, 15)) | id,
           na.action = na.exclude, control = lmeControl(opt = "optim"))

fm3 <- mixed_model(ascites ~ year * sex, data = pbc2,
                   random = ~ year || id, family = binomial())
Mixed <- list(fm1, fm2, fm3)
Cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

#system.time(obj <- jm(Cox, Mixed, time_var = "year"))

#summary(obj)
#traceplot(obj, "alphas")
#gelman_diag(obj)

fForms <- ~ value(log(serBilir)):I(1 * (sex == 'male')) + value(prothrombin) + value(ascites)
fForms <- ~ value(log(serBilir)) + vsquare(value(log(serBilir))) +
    vcubic(value(log(serBilir))) + slope(log(serBilir)) +
    value(log(serBilir)):slope(log(serBilir))
fForms <- ~ value(log(serBilir)) * slope(prothrombin)

# system.time(obj <- jm(Cox, Mixed, time_var = "year", functional_forms = fForms))


Surv_object = Cox
Mixed_objects = Mixed
time_var = 'year'
functional_forms = NULL#fForms
which_independent = cbind(1, 2)
recurrent = FALSE
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con
control$n_chains = 1

##########################################################################################
##########################################################################################

# Cox model for the composite event death or transplantation
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)

# a linear mixed model for log serum cholesterol with random ordering
set.seed(123)
pbc2_ro <- pbc2[sample(seq_len(nrow(pbc2)), nrow(pbc2)), ]
fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
fm2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)

jointFit1 <- jm(CoxFit, list(fm1), time_var = "year",
                functional_forms = ~ area(log(serBilir)))
jointFit2 <- jm(CoxFit, list(fm1), time_var = "year",
                functional_forms = ~ area(log(serBilir), time_window = 5))

summary(jointFit1)
summary(jointFit2)



Surv_object = CoxFit
Mixed_objects = list(fm1, fm2)
time_var = 'year'
functional_forms = ~ area(log(serBilir)) + area(prothrombin)
recurrent = FALSE
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#

time = st
terms = terms_FE_noResp
data = dataL
timeVar = time_var
xxx1 <- degn_matr_area(time, terms, Xbar, time_window)
xxx2 <- degn_matr_area2(time, terms, Xbar, time_window)

all.equal(xxx1, xxx2)

degn_matr_area2 <- function (time, terms, Xbar, time_window) {
    if (!is.list(time)) {
        time <- if (is.matrix(time)) split(time, row(time))
        else split(time, seq_along(time))
    }
    GK <- gaussKronrod(15L)
    wk <- GK$wk
    sk <- GK$sk
    quadrature_points <- function (x, time_window) {
        if (is.null(time_window)) {
            P <- unname(x / 2)
            sk <- outer(P, sk + 1)
            # we divide with x to obtain the area up to time t, divided by t
            # to account for the length of the interval
            list(P = c(t(outer(P / x, wk))), sk = sk)
        } else {
            P <- unname(c(x - x + time_window) / 2)
            sk <- outer(P, sk) + (c(x + x - time_window) / 2)
            # we divide with (x - time_window) to obtain the area from time_window
            # up to time t, divided by t - time_window to account for the length
            # of the interval
            list(P = c(t(outer(P / time_window, wk))), sk = sk)
        }
    }
    sum_qp <- function (m) {
        n <- nrow(m)
        grp <- rep(seq_len(round(n / 15)), each = 15L)
        rowsum(m, grp, reorder = FALSE)
    }
    K <- length(terms)
    out <- vector("list", K)
    for (i in seq_len(K)) {
        time_window_i <- if (length(time_window[[i]])) time_window[[i]][[1L]] else NULL
        terms_i <- terms[[i]]
        qp <- lapply(time, quadrature_points, time_window = time_window_i)
        ss <- lapply(qp, function (x) c(t(x[['sk']])))
        Pwk <- unlist(lapply(qp, '[[', 'P'), use.names = FALSE)


        D <- LongData_HazardModel(ss, data, data[[timeVar]],
                                  data[[idVar]], timeVar,
                                  match(idT, unique(idT)))
        mf <- model.frame.default(terms_i, data = D)
        X <- Pwk * model.matrix.default(terms_i, mf)
        out[[i]] <- sum_qp(X)
    }
    out
}


##########################################################################################
##########################################################################################

fm1 <- lme(log(serBilir) ~ year * sex,
           data = pbc2, random = ~ year | id)
fm2 <- lme(serChol ~ ns(year, 3) + sex + age, data = pbc2, random = ~ year | id,
           na.action = na.exclude)
fm3 <- mixed_model(hepatomegaly ~ sex + age + year, data = pbc2,
                   random = ~ year | id, family = binomial())
fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ year | id, family = binomial())

CoxFit <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

#CoxFit <- survreg(Surv(years, yearsU, status3, type = "interval") ~ 1,
#                  data = pbc2.id, model = TRUE)

fForms <- list("log(serBilir)" = ~ value(log(serBilir)) + slope(log(serBilir)) +
                   slope(log(serBilir)):sex,
               "serChol" = ~ value(serChol) + slope(serChol),
               "hepatomegaly" = ~ vexpit(value(hepatomegaly)) + sex,
               "ascites" = ~ Dexpit(value(ascites)):slope(ascites) + area(ascites))

test <- jm(CoxFit, list(fm1, fm2, fm3, fm4), time_var = "year",
           functional_forms = fForms)

tt <- jm(CoxFit, fm1, time_var = "year", parallel = "multicore")

####

Surv_object = CoxFit
Mixed_objects = list(fm1, fm2, fm3, fm4)
time_var = "year"
functional_forms = fForms
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL


################################################################################
################################################################################


pbc <- pbc2[c("id", "serBilir", "drug", "year", "years", "status2", "spiders")]
pbc$start <- pbc$year
splitID <- split(pbc[c("start", "years")], pbc$id)
pbc$stop <- unlist(lapply(splitID,
                          function(d) c(d$start[-1], d$years[1]) ))

pbc$event <- with(pbc, ave(status2, id,
                           FUN = function(x) c(rep(0, length(x) - 1), x[1])))
pbc <- pbc[!is.na(pbc$spiders), ]
pbc <- pbc[pbc$start != 0, ]


fm1 <- lme(log(serBilir) ~ drug * ns(year, 2), data = pbc, random = ~ ns(year, 2) | id)

tdCox.pbc <- coxph(Surv(start, stop, event) ~ drug * spiders, data = pbc)


####

Surv_object = tdCox.pbc
Mixed_objects = fm1
time_var = 'year'
functional_forms = NULL
recurrent = FALSE
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con
control$n_chains = 1

################################################################################

fm <- lme(sqrt(CD4) ~ obstime * drug, data = aids, random = ~ obstime | patient)
gm <- coxph(Surv(Time, death) ~ drug, data = aids.id)
jmFit <- jm(gm, fm, time_var = "obstime",
            functional_forms = ~ value(sqrt(CD4)) + square(value(sqrt(CD4))))
summary(jmFit)

traceplot(jmFit)

Surv_object = gm
Mixed_objects = fm
time_var = "obstime"
functional_forms = ~ value(sqrt(CD4)) + square(value(sqrt(CD4)))
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con




################################################################################

pbc2$serCholD <- as.numeric(pbc2$serChol > 250)
fm <- mixed_model(serCholD ~ year + age, data = pbc2,
                  random = ~ year | id, family = binomial())
Cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

jFit1 <- jm(Cox, fm, time_var = "year", n_iter = 6500L, n_burnin = 2000L,
            functional_forms = ~ value(serCholD))

jFit2 <- jm(Cox, fm, time_var = "year", n_iter = 6500L, n_burnin = 2000L,
            functional_forms = ~ vexpit(value(serCholD)))

jFit3 <- jm(Cox, fm, time_var = "year", n_iter = 6500L, n_burnin = 2000L,
            functional_forms = ~ vexpit(value(serCholD)) + Dexpit(slope(serCholD)))

summary(jFit1)
summary(jFit2)
summary(jFit3)

traceplot(jFit1, "alphas")
traceplot(jFit2, "alphas")
traceplot(jFit3, "alphas")


fForms <- list("ascites" = ~ vexpit(value(ascites)))
fForms <- list("ascites" = ~ vexpit(value(ascites)):slope(ascites) +
                   value(ascites)*sex + area(ascites) + slope(ascites))

pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
fm1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
fm2 <- lme(prothrombin ~ year * sex, data = pbc2, random = ~ year | id)
fm3 <- mixed_model(ascites ~ year + sex, data = pbc2,
                   random = ~ year | id, family = binomial())

fForms <- ~ log(value(log(serBilir)))

jointFit3 <- update(jointFit2, functional_forms = fForms)
summary(jointFit3)


Surv_object = CoxFit
Mixed_objects = list(fm1, fm2, fm3)
time_var = 'year'
functional_forms = fForms
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con
control$n_chains = 1


fm <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
jm1 <- jm(CoxFit, fm, time_var = "year", cores = 1L)
jm2 <- jm(CoxFit, fm, time_var = "year", cores = 2L)
jm3 <- jm(CoxFit, fm, time_var = "year", cores = 3L)

summary(jm1)
summary(jm2)
summary(jm3)

length(jm1$mcmc$bs_gammas)
length(jm2$mcmc$bs_gammas)
length(jm3$mcmc$bs_gammas)


##############################################################################

data("pbc2", "pbc2.id", package = "JMbayes2")
pbc2$id_char <- as.character(pbc2$id)
pbc2.id$id_char <- as.character(pbc2.id$id)

pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
fm <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id_char)

Surv_object = CoxFit
Mixed_objects = fm
time_var = 'year'
functional_forms = NULL
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
recurrent = FALSE
#



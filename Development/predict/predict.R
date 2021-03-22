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
    source("./R/predict_funs.R")
    Rcpp::sourceCpp('src/mcmc_fit.cpp')

}


object <- jointFit1
ND <- pbc2[pbc2$id %in% c(2, 3, 15), ]
ND$id <- factor(ND$id)
ND2 <- ND[ND$year > 1, ]
ND <- ND[ND$year < 1, ]
ND$status2 <- 0
ND$years <- with(ND, ave(year, id, FUN = function (x) max(x, na.rm = T)))
ND$prothrombin[c(3, 5, 8)] <- NA
newdata = ND
newdata2 = ND2
process = "event"
type_pred = "response"
type = "subject_specific"
level = 0.95; return_newdata = FALSE
n_samples = 200L; n_mcmc = 25L; cores = NULL; seed = 123L

#############################################################
#############################################################

if (is.null(cores)) {
    n <- length(unique(newdata[[object$model_info$var_names$idVar]]))
    cores <- if (n > 20) 4L else 1L
}
components_newdata <- get_components_newdata(object, newdata, n_samples,
                                             n_mcmc, cores, seed)







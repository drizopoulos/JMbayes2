library("survival")
library("nlme")
library("GLMMadaptive")
library("coda")
library("lattice")
library("splines")
source("./R/jm.R")
source("./R/jm_fit.R")
source("./R/help_functions.R")
source("./R/basic_methods.R")
source("./R/create_Wlong_mats.R")
source("./Development/jm/PBC_data.R")

fm1 <- lme(log(serBilir) ~ year * drug + sex + I(year^2) + age + sex:year,
           data = pbc2, random = ~ year | id)
fm2 <- lme(prothrombin ~ year + sex, data = pbc2, random = ~ year | id)
pbc2_nas <- pbc2
pbc2_nas$prothrombin[pbc2_nas$id %in% paste0("A", 1:5)] <- NA
gm1 <- lme(log(serBilir) ~ year * drug + sex + I(year^2) + age + sex:year,
           data = pbc2_nas, random = ~ year | id)
gm2 <- lme(prothrombin ~ year + sex, data = pbc2_nas, random = ~ year | id,
           na.action = na.exclude)
Cox <- coxph(Surv(years, status2) ~ 1, data = pbc2.id)

###############################################################################
###############################################################################

Mixed1 <- list(fm1, fm2)
Mixed2 <- list(gm1, gm2)

Surv_object = Cox
Mixed_objects = Mixed2
time_var = "year"
functional_forms = NULL
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con
control$n_chains = 1
model_info$family_names <- sapply(model_info$families, "[[", "family")
model_info$links <- sapply(model_info$families, "[[", "link")
model_data$y[] <- lapply(model_data$y, as.matrix)
binomial_data <- model_info$family_names == "binomial"
trials_fun <- function (y) {
    if (NCOL(y) == 2) y[, 2] <- y[, 1] + y[, 2]
    y
}
model_data$y[binomial_data] <-
    lapply(model_data$y[binomial_data], trials_fun)

###############################################################################
###############################################################################

lL <- function (model_data, model_info, initial_values) {
    y <- model_data[["y"]]
    X <- model_data[["X"]]
    Z <- model_data[["Z"]]
    betas <- initial_values[["betas"]]
    b <- initial_values[["b"]]
    sigmas <- exp(initial_values[["log_sigmas"]])
    idL <- model_data[["idL"]]
    idL_lp <- model_data[["idL_lp"]]
    unq_idL <- model_data[["unq_idL"]]
    ##
    K <- length(y)
    n <- max(sapply(unq_idL, length))
    out <- numeric(n)
    for (i in seq_len(K)) {
        eta_i <- c(X[[i]] %*% betas[[i]]) +
            rowSums(Z[[i]] * b[[i]][idL_lp[[i]], , drop = FALSE])
        log_y <- dnorm(y[[i]], eta_i, sigmas[[i]], log = TRUE)
        out[unq_idL[[i]]] <- out[unq_idL[[i]]] +
            rowsum(log_y, idL_lp[[i]], reorder = FALSE)
    }
    list(log_Lik_vec = out, log_Lik = sum(out))
}

logLong_R <- lL(model_data, model_info, initial_values)
logLong_cpp <- test_log_long(model_data, model_info, initial_values)

logLong_R$log_Lik
logLong_cpp$log_Lik

all.equal(logLong_R$log_Lik_vec, c(logLong_cpp$log_Lik_vec))



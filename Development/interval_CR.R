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
Rcpp::sourceCpp('src/mcmc_fit.cpp')
load("C:/Users/drizo/OneDrive/Desktop/JMbayes2 test/trainingset_1.RData")
load("C:/Users/drizo/OneDrive/Desktop/JMbayes2 test/trainingset_id_1.RData")
load("C:/Users/drizo/OneDrive/Desktop/JMbayes2 test/Joint model_1.RData")
################################################################################
################################################################################

f <- function (x) {
    delta <- x$status.cmp
    if (delta == 0) {
        out <- x[c(1, 1), ]
        out$start <- 0
        out$stop <- c(x$time.cmp2, x$time.cmp1)
        out$event <- c(0, 0)
        out$weight <- 1
        out$strt <- c("trt", "prg")
        out$intgr_ <- FALSE
    } else if (delta == 2) {
        out <- x[c(1, 1), ]
        out$start <- 0
        out$stop <- c(x$time.cmp2, x$time.cmp1)
        out$event <- c(1, 0)
        out$weight <- 1
        out$strt <- c("trt", "prg")
        out$intgr_ <- FALSE
    } else if (delta == 1) {
        out <- x[rep(1, 5), ]
        out$start <- 0
        a <- x$time.cmp1
        b <- x$time.cmp2
        out$stop <- c(x$time.cmp2, a, (2 * a + b) / 3, (a + 2 * b) / 3, b)
        n1 <- nrow(out) - 1
        out$event <- c(0, rep(1, n1))
        out$weight <- c(1, rep((b - a) / 8, n1))
        out$strt <- c("trt", rep("prg", n1))
        out$intgr_ <- c(FALSE, rep(TRUE, n1))
    }
    out
}
eventDF <- do.call("rbind", lapply(split(train.dat.id, train.dat.id$CISNET_ID), f))
row.names(eventDF) <- seq_len(nrow(eventDF))
eventDF$strt <- factor(eventDF$strt, levels = c("trt", "prg"))
eventDF$stop[eventDF$start == eventDF$stop] <- 1e-05
attr(eventDF$weight, "integrate") <- eventDF$intgr_


CoxFit <- coxph(Surv(start, stop, event) ~ density:strata(strt), data = eventDF,
                weights = eventDF$weight, model = TRUE)

lmeFit <- lme(PSAValue ~ ns(time, k = c(1.49, 3.535), B = c(0, 10.223)) + DxAge,
              data = train.dat,
              random = list(CISNET_ID = pdDiag(form = ~ ns(time, k = c(1.49, 3.535), B = c(0, 10.223)))))

lmeFit <- lme(PSAValue ~ ns(time, df = 3) + DxAge,
              data = train.dat, control = lmeControl(opt = 'optim'),
              random = ~ ns(time, df = 3) | CISNET_ID)

jmFit <- jm(CoxFit, lmeFit, time_var = "time",
            functional_forms = ~ value(PSAValue):strt +
                slope(PSAValue,eps = 1, direction = "back"):strt,
            n_iter = 6500L, n_burnin = 500L)
summary(jmFit)

iccsjm.model$summary$coef



Surv_object = CoxFit
Mixed_objects = lmeFit
time_var = 'time'
functional_forms = NULL#fForms
which_independent = NULL
recurrent = FALSE
data_Surv = NULL
id_var = NULL
priors = NULL
control = NULL
#
model_data <- Data
control <- con
control$n_chains = 1


xx <- rnorm(length(model_data$id_h2))
lse <- function (xx, id_h2, intgr_ind) {
    unq_idh2 <- unique(id_h2)
    out <- numeric(length(unq_idh2))
    for (i in unq_idh2) {
        idx <- id_h2 == i
        intgr_i <- intgr_ind[idx]
        nn <- length(intgr_i)
        if (nn > 0) {
            out[i] <- logSumExp(xx[idx])
        } else {
            out[i] <- xx[idx]
        }
    }
    out
}

all.equal(lse(xx, id_h2, intgr_ind), c(lseC(xx, id_h2 - 1, intgr_ind)))

str(lseC(xx, model_data$id_h2 - 1, model_data$intgr_ind))


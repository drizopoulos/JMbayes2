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
source("./R/create_Wlong_mats.R")
load("./Development/CI/model.RData")
source("./Development/CI/prepare_data.R")
Rcpp::sourceCpp('src/mcmc_fit.cpp')


newdataL <- pbc2[pbc2$id == 14, ]
newdataE <- pbc2_CR[pbc2_CR$id == 14, ]
newdataE$event <- 0
newdata <- list(newdataL = newdataL, newdataE = newdataE)

object = jointFit
newdata = newdata
newdata2 = NULL
times = NULL
process = "event"
type_pred = "response"
type = "subject_specific"
level = 0.95; return_newdata = TRUE
n_samples = 200L; n_mcmc = 55L; cores = NULL
seed = 123L


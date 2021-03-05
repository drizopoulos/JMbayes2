# packages
library(JMbayes)
library(splines)
library(mstate)
require(MASS)
require(Matrix)

# Source sim function
source(file.path(getwd(), "Development/Dev_Local_GP/MS_CR/sim_ms_basic.R"))

simdata <- sim_mstate(seed = 2021)

save(simdata, file = file.path(getwd(), "Development/Dev_Local_GP/MS_CR/simdata_ms_vignette.RData"))

# Simulate 100 datasets
#simdats <- list(NULL)
#for(i in 1:50) {
#  simdats[[i]] <- sim_mstate(seed = i)
#  print(paste0("dataset ", i, " done"))
#}

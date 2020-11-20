library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(ggplot2)

#load(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\jm_manual_2808.Rdata')
load(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/passed_to_mcmc_fit.RData')

sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\log_u.cpp')
sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\log_long.cpp')
sourceCpp(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Dev_Local\\Cpp_funs\\Xbetas_calc.cpp')

Xbetas <- Xbeta_calc(model_data$X, initial_values$betas)
scales <- rep(0.1, 312)
extra_parms <- rep(0, 4)
families_chr <- c('gaussian', 'gaussian', 'binomial', 'binomial')
links_chr <- c('identity', 'identity', 'logit', 'logit')

idL_lp[[1]] <- idL_lp[[1]] - 1
idL_lp[[2]] <- idL_lp[[2]] - 1
idL_lp[[3]] <- idL_lp[[3]] - 1
idL_lp[[4]] <- idL_lp[[4]] - 1

idllp <- idL_lp

idllp[[1]] <- c(311)

eta <- rep(1, 312)

JM1$model_data$y[] <- lapply(JM1$model_data$y, as.matrix)
model_data$idL

log_long(model_data$y, eta, scales, extra_parms, families_chr, links_chr, 
         model_data$idL, model_data$unq_idL)

log_u(Xbetas, model_data$Z, initial_values$b, lapply(model_data$idL_lp, FUN = function(x) x-1 ), 
      model_data$y, scales, extra_parms, model_info$family_names, model_info$links, 
      lapply(model_data$idL, FUN = function(x) x-1 ), lapply(model_data$unq_idL, FUN = function(x) x-1 ))


b[[3]][idL_lp[[3]]]

lapply(idL_lp, range)

class(JM1$model_data$y[[1]])

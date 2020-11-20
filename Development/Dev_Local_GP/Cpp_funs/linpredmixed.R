library(Rcpp)
library(RcppArmadillo)

sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/linpredmixed.cpp')

obj$model_data$idL_lp_minus1 <- lapply(obj$model_data$idL_lp, FUN = function(x) x - 1)
eta <- linpred_mixed(obj$model_data$X, obj$initial_values$betas, obj$model_data$Z, obj$initial_values$b, obj$model_data$idL_lp_minus1)

eta[[1]]
length(obj$model_data$idL[[1]])

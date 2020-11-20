library(Rcpp)
library(RcppArmadillo)

sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Cpp_funs/idL_lp_fast.cpp')

chk <- foo(obj$model_data$idL)
chk[, 1]

split(obj$model_data$idL, obj$model_data$idL)
lapply(obj$model_data$idL, FUN = function(x) tapply(x, x, length))

klain1 <- seq_along(obj$model_data$idL)

klain1 <- lapply(obj$model_data$idL, function(x)  seq_along(x) - 1)
klain2 <- mapply(function(x, y) split(x, y), klain1, obj$model_data$idL, SIMPLIFY = TRUE)

klain2[[1]][[1]]

eta_new = lapply(eta, function (x) x * 0)

klain3 <- foo(eta, eta_new, klain2)


eta[[1]][klain2[[1]]]

klain3[1, ]

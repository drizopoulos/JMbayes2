save(res_thetas, statistics, model_data, model_info, control, 
     file = paste0(getwd(), '/Development/Dev_Local_GP/update/update_new/input_mloglik_new_default.RData'))

out$`1`$mcmc$cumsum_b

library(Rcpp)     
library(RcppArmadillo)

Rcpp::sourceCpp(file = paste0(getwd(), '/Development/Dev_Local_GP/update/update_new/chck_burnin.cpp'))

check(0, 0)
0 + 1 + 2 + 3 + 4 + 5 + 6 + 7 + 8 + 9
dim(out$`1`$mcmc$b)
klain1 <- out$`1`$mcmc$b[,,1]
for(i in 2:3000) {
  klain1 <- klain1 + out$`1`$mcmc$b[,,i]
}

all.equal(klain1, out$`1`$mcmc$cumsum_b)


out$`1`$mcmc$cumsum_b
dim(out$`1`$mcmc$bs_gammas)

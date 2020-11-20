# package
library(Rcpp)
library(RcppArmadillo)

# load cpp environment
load(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Eta_with_NAs/passed_to_cpp.RData')

# source cpp
sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Eta_with_NAs/create_eta.cpp')
sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Eta_with_NAs/eta_check.cpp')
#sourceCpp(file = 'D:/gpapageorgiou/Github_Repos/JMbayes2/Dev_Local/Eta_with_NAs/field_3d.cpp')

# create example eta
eta <- linpred_mixed(model_data, model_info, initial_values, priors, control, vcov_prop)
eta_prop <- eta
for (i in 1:4) {
  for(j in 1:length(eta_prop[[i]])) {
    eta_prop[[i]][[j]] <- eta_prop[[i]][[j]] - eta_prop[[i]][[j]]
  }
}

unq_idL_lst <- lapply(model_data$unq_idL, function(x) split(x, x))
unq_idL_check <- unique(unlist(lapply(unq_idL_lst, names)))
unq_idL_lst <- setNames(unq_idL_lst, 1:length(unq_idL_lst))
unq_idL_lst <- lapply(do.call(mapply, c(FUN = c, lapply(unq_idL_lst, `[`, unq_idL_check))), function(x) as.numeric(names(x)))
unq_idL_lst[[41]]

# check stack overflow
# https://stackoverflow.com/questions/18538977/combine-merge-lists-by-elements-names

eta[[1]]
eta_prop[[1]]
length(41:312)
#chk <- sample(c(0, 1), 312, replace = T)
chk <- c(rep(0, 39), 1, 1, rep(0, 271))
klain <- lapply(lapply(model_data$idL, function(x) x - 1), function(x) split(x, x))
lapply(klain, function(x) unique(names(x)))
klain[[2]]

chk1 <- eta_check_wrap(eta, eta_prop, chk, klain, 312, model_data, unq_idL_lst)

chk1[[2]][[12]]
eta[[1]][[3]]
model_data$idL[[1]] == unique(klain[[1]][[2]]) + 1

model_data$idL[[1]] == 2

chk_indx_39_1 <- which(model_data$idL[[1]] == 39)
chk_indx_39_2 <-  which(model_data$idL[[2]] == 39)
chk_indx_40_1 <- which(model_data$idL[[1]] == 40)
chk_indx_40_2 <-  which(model_data$idL[[2]] == 40)
chk_indx_41_1 <- which(model_data$idL[[1]] == 41)
chk_indx_41_2 <-  which(model_data$idL[[2]] == 41)
chk_indx_42_1 <- which(model_data$idL[[1]] == 42)
chk_indx_42_2 <-  which(model_data$idL[[2]] == 42)



chk1[[1]][chk_indx_39_1]
chk1[[1]][chk_indx_40_1]
chk1[[1]][chk_indx_41_1]
chk1[[1]][chk_indx_42_1]

chk1[[2]][chk_indx_39_2]
chk1[[2]][chk_indx_40_2]
chk1[[2]][chk_indx_41_2]
chk1[[2]][chk_indx_42_2]





####################

# create List of Lists of uvec version of idL
# directly with c++ indexing
model_data$idL_LstOfLst <- lapply(lapply(model_data$idL, function(x) x - 1), function(x) split(x, x))
# create a new list which
# has elements equal to the number of subjects and 
# each element represents a subject
# each element is a uvec
# each uvec indicates for which longitudinal outcomes 
# does the corresponding subjects have observations for
unq_idL_outc_lst <- lapply(model_data$unq_idL, function(x) split(x, x))
unq_idL_outc_names <- unique(unlist(lapply(unq_idL_outc_lst, names)))
unq_idL_outc_lst <- setNames(unq_idL_outc_lst, 1:length(unq_idL_outc_lst))
unq_idL_outc_lst <- lapply(do.call(mapply, c(FUN = c, lapply(unq_idL_outc_lst, `[`, unq_idL_outc_names))), 
                           function(x) as.numeric(names(x)))
model_data$unq_idL_outc_lst <- unq_idL_outc_lst

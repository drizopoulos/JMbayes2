load(file = 'D:\\gpapageorgiou\\Github_Repos\\JMbayes2\\Development\\Dev_Local_GP\\update\\update_wrkng\\example_objects_inits.RData')

names(initial_values)

object1 <- joint_model_fit_1
object2 <- joint_model_fit_2
object3 <- joint_model_fit_3
object4 <- joint_model_fit_4

dim(joint_model_fit_3$mcmc$b[[1]])

# bs_gammas:        6000 12
# tau_bs_gammas:    6000  1
# gammas:           6000  1
# alphas:           6000  1
# W_std_gammas:     6000  1
# Wlong_std_alphas: 6000  1
# D                 6000  3
# betas1:           6000  7
# sigmas1:          6000  1
# b:                312   2  1 OR 312 2 6000

klain <- extract_last_iteration(object1$mcmc)
klain$bs_gammas$`1` $bs_gammas[[1]]

klain$bs_gammas$`2`
dim(klain$b[[1]])
initial_values$bs_gammas

extract_initial_values <- function(object) {
  if(!inherits(object, 'jm'))
    stop("Use only with 'jm' object. \n")
  #out <- list(NULL)
  #dimb <- unique(do.call(c, lapply(object$mcmc$b, function(x) dim(x)[length(dim(x))])))
  #if (dimb == 1) {
  #  joint_model_fit_1$mcmc$betas1
  #} else {
    
  #}
  
}

extract_initial_values <- function(object) {
  if(!inherits(object, 'jm'))
    stop("Use only with 'jm' object. \n")
  last_iters <- extract_last_iteration(object$mcmc)
  
}

dim(object1$mcmc$b$`1`)[length(dim(object1$mcmc$b$`1`), n = 1L)]


extract_init_vals <- function(object) {
  if (!inherits(object, 'jm'))
    stop("Use only with 'jm' object.\n")
  #initial_values <- object$statistics$Mean
  
}
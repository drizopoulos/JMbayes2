#rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
setwd("./../..")

source(file.path(getwd(), "Development/MCMC/Surv_Model/sample_Surv_Funs.R")) # log_density_surv(), logPrior(), rmvnorm()
source(file.path(getwd(), "Development/MCMC/Surv_Model/examples.R")) 

target_log_dist <- function (betas, i) {
  
  linear_predictor <- linpred_mixed(X = X[i], 
                                    betas = betas[i], 
                                    Z = Z[i], 
                                    b = b[i], 
                                    id = idL_lp[i]) 
  
  log_Y <- log_density_mixed(y = y[i],
                             linear_predictor = linear_predictor, # No need for [i], because the linear_predictor is a list of length 1. We don't need the likelihood from all mixed-models here. The posterior is proportional to its own likelihood
                             log_sigmas = log_sigmas[i], 
                             Funs = Funs[i],
                             mu_funs = mu_funs[i],
                             nY = n,
                             unq_idL = unq_idL[i],
                             idL_lp = idL[i]) #? or should it be idL_lp[i]? # > sum(!idL_lp[[2]]==idL[[2]]) [1] 975

  log_S <- log_density_surv(bs_gammas = bs_gammas,
                            gammas = gammas,
                            alphas = alphas) 
  
  log_prior <- logPrior(theta = betas[[i]],
                        mean_theta = prior_mean_betas[[i]],
                        Tau_theta = prior_Tau_betas[[i]],
                        tau_theta = 1)
  
  sum(log_Y) + log_S + log_prior  #? do I need to do the sum(log_Y), or am I doing something wrong?
  
}

# Parameters -------------------------------------------------------------------

Funs <- lapply(test$model_info$families, log_dens_Funs)
mu_funs <- lapply(test$model_info$families, "[[", 'linkinv')
n <- test$model_data$n
unq_idL <- test$model_data$unq_idL
idL_lp <- test$model_data$idL_lp
idL <- test$model_data$idL

# parameter values
betas <- test$initial_values$betas
b <- test$initial_values$b
gammas <- test$initial_values$gammas
bs_gammas <- test$initial_values$bs_gammas
alphas <- test$initial_values$alphas
log_sigmas <- test$initial_values$log_sigmas

# outcome vectors and design matrices
X <- test$model_data$X
Z <- test$model_data$Z
y <- test$model_data$y
W0_H <- test$model_data$W0_H
W_H <- test$model_data$W_H
#? the code below seems to be missing from examples.R, I copy-pasted it from prepare_test_fast_log_surv.R
# id_H is used to repeat the random effects of each subject GK_k times
id_H <- lapply(X_H, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
# this is the linear predictor for the longitudinal outcomes evaluated at the
# Gauss-Kronrod quadrature points
eta_H <- linpred_surv(X_H, betas, Z_H, b, id_H)
# Wlong is the design matrix of all longitudinal outcomes according to the specified
# functional forms per outcome already multiplied with the interaction terms matrix U
Wlong_H <- create_Wlong(eta_H, functional_forms_per_outcome, U_H)
if (length(which_event)) {
  id_h <- lapply(X_h, function (x) seq_len(nrow(x[[1]])))
  eta_h <- linpred_surv(X_h, betas, Z_h, b, id_h)
  Wlong_h <- create_Wlong(eta_h, functional_forms_per_outcome, U_h)
}
if (length(which_interval)) {
  id_H2 <- lapply(X_H2, function (i, n) rep(seq_len(n), each = control$GK_k), n = n)
  eta_H2 <- linpred_surv(X_H2, betas, Z_H, b, id_H2)
  Wlong_H2 <- create_Wlong(eta_H2, functional_forms_per_outcome, U_H2)
} else {
  Wlong_H2 <- rep(list(matrix(0.0, length(Time_right), 1)), length(W_H))
}


# priors
prior_mean_betas <- test$priors$mean_betas
prior_Tau_betas <- test$priors$Tau_betas

vcov_prop_betas <- lapply(list(fm1, fm2, fm3, fm4), vcov2)
#vcov_prop_betas <- test$vcov_prop$vcov_prop_betas
#vcov_prop_tilde_betas <- test$vcov_prop$vcov_prop_tilde_betas

# Sampling ---------------------------------------------------------------------

M <- 5000L # number of iterations
nburnin <- 500L 
#? thinning parameter?

current_betas <- betas
acceptance_betas <- matrix(0.0, M, length(betas)) # No need for a list/array because we're block-updating for each mixed-model
res_betas <- lapply(betas, function(b){ m <- matrix(0.0, M, length(b))
                                        colnames(m) <- names(b)
                                        m})

scale_betas <- rep(0.1, length(betas))

# Update betas
system.time({#set.seed(2020)
for (i in seq_along(betas)) { # i-th mixed model

  for (m in seq_len(M)) { # m-th sample
    
    if (m == 1) denominator <- target_log_dist(current_betas, i) # + 0 given the proposal is symmetric
    
      proposed_betas <- current_betas
      
      proposed_betas[[i]] <- rmvnorm(n = 1,
                                     mu = current_betas[[i]],
                                     Sigma = vcov_prop_betas[[i]] * scale_betas[i])
      
      numerator <- target_log_dist(proposed_betas, i) # + 0 given the proposal is symmetric
      
      log_ratio <- numerator - denominator
      
      if (log_ratio > log(runif(1))) {
        current_betas <- proposed_betas
        acceptance_betas[m, i] <- 1.0
        denominator <- numerator
      }
      
      if (m > 20) {
        scale_betas[i] <- robbins_monro_univ(scale = scale_betas[i],
                                             acceptance_it = acceptance_betas[m, i],
                                             it = m, target_acceptance = 0.234)
      }
      
      res_betas[[i]][m, ] <- current_betas[[i]]
      
    }

}
})

#user  system elapsed 
#13.05    0.05   13.11

# Results ----------------------------------------------------------------------

formulas <- c("lme(log(serBilir) ~ year * sex + I(year^2) + age + prothrombin, data = pbc2, random = ~ year | id)",
              "lme(serChol ~ ns(year, 3) + sex + age, data = pbc2, random = ~ year | id, na.action = na.exclude)",
              "mixed_model(hepatomegaly ~ sex + age, data = pbc2, random = ~ 1 | id, family = binomial())",
              "mixed_model(ascites ~ year + age, data = pbc2, random = ~ 1 | id, family = binomial())")

for (i in seq_along(current_betas)) { # i-th mixed model
  
  samples <- res_betas[[i]]
  ncol <- 2
  nrow <- ceiling(ncol(samples)/2)
  par(mfrow = c(nrow, ncol))
  
  for(j in seq_len(ncol(samples))) { # j-th parameter
    plot(seq_len(M), samples[, j], type = 'l',
         xlab= "Iterations", ylab= "Samples", bty = "l")
    lines(seq_len(nburnin),samples[seq_len(nburnin), j], col = "darkgrey")
    abline(h= mean(samples[-seq_len(nburnin), j]), 
           lty= 2, col= "red") # post mean
    abline(h= quantile(samples[-seq_len(nburnin), j], probs = c(0.025, 0.975)), 
           lty= 3, col= "red") # post 95CI
    title(colnames(samples)[j], line = 0)
  }
  
  ar <- round(mean(acceptance_betas[-seq_len(nburnin), i]), 2)
  mtext(text = paste0("fm", i, " (acceptance rate ", ar, ")"), 
        side = 3, line = -1.5, outer = TRUE, font = 2)
  mtext(text = formulas[i], side = 3, line = -3, outer = TRUE)
  
}

# optimal accep. rates MULTIvariate proposal: 0.234 (Roberts et al., 1997)
# optimal accep. rates UNIvariate proposal: 0.44 (Roberts and Rosenthal, 2001)

# fm1 <- lme(log(serBilir) ~ year * sex + I(year^2) + age + prothrombin,
#            data = pbc2, random = ~ year | id)
# fm2 <- lme(serChol ~ ns(year, 3) + sex + age, data = pbc2, random = ~ year | id,
#            na.action = na.exclude)
# fm3 <- mixed_model(hepatomegaly ~ sex + age, data = pbc2,
#                    random = ~ 1 | id, family = binomial())
# fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
#                    random = ~ 1 | id, family = binomial())
# 
# CoxFit <- coxph(Surv(years, status2) ~ age * sex,
#                 data = pbc2.id, model = TRUE)
klain <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)

klain[length(klain)]

lapply(initial_values, class)
lapply(klain, class)

names(klain)

unlist(klain$bs_gammas)

?parallel::parLapply
klain
initial_values$betas

initial_values$log_sigmas


# data 
data(Boston)

# function - calculate the mse from a model fit on bootstrapped samples from the Boston dataset
model.mse <- function(x) {
  id <- sample(1:nrow(Boston), 200, replace = T)
  mod <- lm(medv ~ ., data = Boston[id,])
  mse <- mean((fitted.values(mod) - Boston$medv[id])^2)
  return(mse)
}

# initialising the list to pass to the apply functions
x.list <- sapply(1:10000, list)

# detect the number of cores
n.cores <- detectCores()
n.cores


x.list <- list(NULL)
ind.list <- list(NULL)
for (i in 1:1000) {
  x.list[[i]] <- rnorm(100)
  ind.list[[i]] <- sample(1:100, 1)
}
library(parallel)
n.cores <- detectCores() - 1
?mapply
?parallel::parLapply
demek1 <- mapply(FUN = function(x, y) x[y], x.list, ind.list)
demek1[1] 
x.list[[1]][ind.list[[1]]]
cl <- parallel::makeCluster(n.cores)
demek2 <- clusterMap(cl, fun = function(x, y) x[y], x.list, ind.list)
demek1[1]
demek2[1]
x.list[[1]][ind.list[[1]]]

demek3 
stopCluster(cl)
  
  
singlecore <- mapply(l.list,  tail, n = 1L)
cl <- parallel::makeCluster(n.cores)
mcore <- parLapply(cl, x.list, tail, n = 1L)
singlecore[[1000]]
mcore[[1000]]

?clusterMap

lapply(seq_len(object1$control$n_chains), function(x) list('x' = 1:10))

mapply(create_HC_X3, ...)

#

x <- X[[1]]
z <- Z[[1]] 
id <- id1
form <- formula(fm1)
data <- dataL

create_HC_X3 <- function(x, z, id, form, data) {
  
  # functions
  check_tv <- function (x, id) {
    !all(sapply(split(x, id),
                function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
  }
  
  # local vars
  cnams_x <- colnames(x)
  cnams_z <- colnames(z)
  n_res <- ncol(z)
  
  X_HC   <- vector("list", length = n_res)
  mat_HC <- matrix(0, nrow= n_res, ncol = ncol(x), 
                   dimnames = list(cnams_z, cnams_x))
  mat_HC[cbind(which(cnams_z %in% cnams_x), which(cnams_x %in% cnams_z))] <- 1 # x_in_z
  
  # baseline (assumes every model has a random intercept)
  x_notin_z <- which(!cnams_x %in% cnams_z)
  baseline <- x_notin_z[!apply(x[, x_notin_z, drop = FALSE], 2L, check_tv, id= id)]
  X_HC[[1]] <- x[!duplicated(id), baseline, drop = FALSE]
  mat_HC[cbind(1, baseline)] <- 2 # baseline
  
  # remaining RE
  if(n_res > 1){
    for(i in seq_len(n_res)[-1]) {

      xint_in_z <- union(grep(paste0(cnams_z[i], ":"), cnams_x), grep(paste0(":", cnams_z[i]), cnams_x)) # interactions can be found as RE:var1, var1:RE, or var1:RE:var2
      if(length(xint_in_z)==0) next
      
      data_temp <- data
      data_temp[[ cnams_z[i] ]] <- 1
      x_temp <- model.matrix(form, data= data_temp)
      
      baseline2 <- xint_in_z[!apply(x_temp[, xint_in_z, drop = FALSE], 2L, check_tv, id= id)]
      X_HC[[i]] <- x_temp[!duplicated(id), baseline2]
      mat_HC[cbind(i, baseline2)] <- 3 # xint_in_z
      
    }
  }
  
  x_in_z_base = which(colSums(mat_HC>0) == 1)

  # return
  list(mat_HC = mat_HC, 
       X_HC = X_HC, 
       x_in_z_base = x_in_z_base,
       nfes_HC = length(x_in_z_base),
       z_in_x = which(rowSums(mat_HC==1) == 1),
       x_in_z = which(colSums(mat_HC==1) == 1),
       x_notin_z = which(colSums(mat_HC) == 0),
       xbas_in_z = mat_HC[, x_in_z_base] > 1,
       baseline = baseline #?? I believe this will not be needed, remove later
       )
  
}

#

#j <- 3
#x <- X[[j]]
#z <- Z[[j]]
#id <- idL[[j]]
#formL <- lapply(Mixed_objects, formula) #?? not sure if we already have this variable in jm()
#form <- formL[[j]]
#data <- dataL

#components_HC <- mapply2(create_HC_X3, X, Z, idL, formL, rep(list(dataL), length(nres)))
#mat_HC <- lapply(components_HC, "[[", "mat_HC")
#X_HC <- lapply(components_HC, "[[", "X_HC")
#x_in_z_base <- lapply(components_HC, "[[", "x_in_z_base")
#nfes_HC <- sapply(components_HC, "[[", "nfes_HC")
#z_in_x <- sapply(components_HC, "[[", "z_in_x")
#x_in_z <- sapply(components_HC, "[[", "x_in_z")
#x_notin_z <- sapply(components_HC, "[[", "x_notin_z")
#xbas_in_z <- sapply(components_HC, "[[", "xbas_in_z")
baseline <- lapply(components_HC, "[[", "baseline")

#X_dot <- create_X_dot3(nres, nfes_HC, z_in_x, x_in_z, X_HC, nT, unq_idL, xbas_in_z)

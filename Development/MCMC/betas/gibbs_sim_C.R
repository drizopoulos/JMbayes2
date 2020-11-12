# Functions
{
 mapply2 <- function (FUN, ..., MoreArgs = NULL, USE.NAMES = TRUE) {
        mapply(FUN, ..., MoreArgs = MoreArgs, SIMPLIFY = FALSE,
               USE.NAMES = USE.NAMES)
      }
      
 .bdiag <- function (mlist) {
      # constructs a block-diagonal matrix
      mlist <- mlist[sapply(mlist, length) > 0]
      #if (length(mlist) == 1)
      #    mlist <- unlist(mlist, recursive = FALSE)
      csdim <- rbind(c(0, 0), apply(sapply(mlist, dim), 1, cumsum))
      ret <- array(0, dim = csdim[length(mlist) + 1, ])
      add1 <- matrix(rep(1:0, 2), ncol = 2)
      for (i in seq_along(mlist)) {
        indx <- apply(csdim[i:(i + 1), ] + add1, 2, function(x) x[1]:x[2])
        if (is.null(dim(indx))) {
          ret[indx[[1]], indx[[2]]] <- mlist[[i]]
        }
        else {
          ret[indx[, 1], indx[, 2]] <- mlist[[i]]
        }
      }
      #colnames(ret) <- unlist(lapply(mlist, colnames))
      ret
    }
  
 create_HC_X2 <- function (x, z, id) {
    check_tv <- function (x, id) {
      !all(sapply(split(x, id),
                  function (z) all(z - z[1L] < .Machine$double.eps^0.5)))
    }
    
    cnams_x <- colnames(x)
    cnams_z <- colnames(z)
    if (!"(Intercept)" %in% cnams_x || !"(Intercept)" %in% cnams_z) {
      stop("cannot perform hierarchical centering in the absense of an ",
           "intercept term in both the fixed and random effects design ",
           "matrices.")
    }
    
    z_in_x <- which(cnams_z %in% cnams_x)
    x_in_z <- which(cnams_x %in% cnams_z)
    x_notin_z <- which(!cnams_x %in% cnams_z)
    baseline <- x_notin_z[!apply(x[, x_notin_z, drop = FALSE], 2L, check_tv, id= id)]
    x_notin_z <- setdiff(x_notin_z, baseline)
    if (!length(baseline)) baseline <- as.integer(NA)
    if (!length(x_notin_z)) x_notin_z <- as.integer(NA)
    list(z_in_x = z_in_x, baseline = baseline, x_in_z = x_in_z, x_notin_z = x_notin_z,
         Xbase = x[!duplicated(id), baseline, drop = FALSE])
  }
  
 create_X_dot <- function(Xbase, ids_unq, ids_out, nres, z_in_x){
  
  # Xbase - a list of the baseline design matrices X per outcome
  # ids_unq - a vector of unique ids
  # ids_out - a list of ids present in each outcome
  # nres - a vector of the number of random effects per outcome
  
  do.call(rbind, lapply(ids_unq, function(id){
    
    Xbase_i <- mapply2(function(out_Xbase, out_ids){
      
      out_Xbase[match(id, unique(out_ids)), ]
      
    }, Xbase, ids_out)
    
    .bdiag(
      mapply2(function(X_outc, diag1_ncol, zinx){
        
        .bdiag(list(t(c(1, X_outc)), diag(diag1_ncol)[, zinx[-1]-1, drop= FALSE]))
        
      }, Xbase_i, nres-1, z_in_x)
    )
    
  }))
  
}
}

# Generata data
{
  set.seed(2020)
  
  n <- 500 # number of subjects
  n_i <- 15 # number of repeated measurements
  tmax <- 10 # maximum follow-up time
  
  data <- data.frame(id    = rep(seq_len(n), each = n_i),
                     time11 = c(replicate(n, c(0, sort(runif(n_i - 1, 1, tmax))))),
                     time12 = c(replicate(n, c(0, sort(runif(n_i - 1, 1, tmax))))),
                     time2 = c(replicate(n, c(0, sort(runif(n_i - 1, 1, tmax))))), 
                     time3 = c(replicate(n, c(0, sort(runif(n_i - 1, 1, tmax))))),
                     group1= rep(rbinom(n, size= 1, prob= 0.5), each= n_i), 
                     group2= rep(rbinom(n, size= 1, prob= 0.5), each= n_i), 
                     group31= rep(rbinom(n, size= 1, prob= 0.5), each= n_i),
                     group32= rep(rbinom(n, size= 1, prob= 0.5), each= n_i)
  )

}

# Add missing data
set.seed(2020)
n_na <- 150 # 150 patients (out of 500) missing per outcome (out of 3)
ids_na <- sample(seq_len(n), size= n_na*3) 
ids_na <- split(ids_na, rep(1:3, times= n_na))

idL <- lapply(ids_na, function(ids_out){ data$id[!data$id %in% ids_out] })

X <- list(model.matrix(~ 1 + time11 + time12 + as.factor(group1), data = data[data$id %in% idL[[1]],]),
          model.matrix(~ 1 + time2 + as.factor(group2), data = data[data$id %in% idL[[2]],]),
          model.matrix(~ 1 + time3 + as.factor(group31) + as.factor(group32) + time3:as.factor(group31), data = data[data$id %in% idL[[3]],]))

Z <- list(model.matrix(~ 1 + time11 + I(time11^2), data = data[data$id %in% idL[[1]],]),
          model.matrix(~ 1 + time2, data = data[data$id %in% idL[[2]],]),
          model.matrix(~ 1, data = data[data$id %in% idL[[3]],]))

nres <- sapply(Z, ncol)
nfes <- sapply(X, ncol)
ind_FE <- split(seq_len(sum(nfes)), rep(seq_along(X), nfes))
ind_RE <- split(seq_len(sum(nres)), rep(seq_along(Z), nres))

componentsHC <- mapply(create_HC_X2, x= X, z= Z, id= idL, SIMPLIFY = FALSE)
x_in_z <- lapply(componentsHC, "[[", "x_in_z")
baseline <- lapply(componentsHC, "[[", "baseline")
x_in_z_base <- mapply2(function (x, y) sort(c(x, y)), x_in_z, baseline)

ind_FE_HC <- unlist(mapply2(function (x, ind) x[ind], ind_FE, x_in_z_base),
                    use.names = FALSE)
 
Xbase <- lapply(componentsHC, "[[", "Xbase")
z_in_x <- lapply(componentsHC, "[[", "z_in_x")
X_dot  <- create_X_dot(Xbase= Xbase, 
                       ids_unq= seq_len(n), #?? change this to account for that fact that some patients might not have longitudinal outcomes
                       ids_out= idL, 
                       nres= nres,
                       z_in_x = z_in_x)



out_in <- sapply(idL, "%in%", x = seq_len(n))
all_pat <- apply(out_in, 1L, paste0, collapse = "/")
id_patt <- match(all_pat, unique(all_pat))
find_patt <- function (patt, n) which(rep(patt, times = n))
nfes_HC <- sapply(x_in_z_base, length)
ind_RE_patt <- apply(unique(out_in), 1L, find_patt, n = nres)
ind_FE_patt <- apply(unique(out_in), 1L, find_patt, n = nfes_HC)

mean_betas_HC <- numeric(sum(nfes_HC)) # priors
Tau_betas_HC <- diag((rep(1/10000, sum(nfes_HC)))) 

betas <- lapply(nfes, numeric)

##

set.seed(2020)
R <- matrix(0, sum(nres), sum(nres))
R[lower.tri(R, diag= FALSE)] <- runif( (sum(nres)^2 - sum(nres))/2, 0.2, 0.8)
R <- R + t(R)
diag(R) <- 1
sd <- runif(sum(nres), 0.5, 3)
D <- diag(sd) %*% R %*% diag(sd)
bb <- MASS::mvrnorm(n, rep(0, sum(nres)), D) 
b <- lapply(ind_RE, function(cols) bb[, cols, drop= FALSE] )


data <- list(n_iter = 100,
             X_dot = X_dot,
             ind_FE_HC = ind_FE_HC,
             id_patt = id_patt,
             ind_RE_patt = ind_RE_patt,
             ind_FE_patt = ind_FE_patt,
             ind_FE = ind_FE,
             mean_betas_HC = mean_betas_HC,
             Tau_betas_HC = Tau_betas_HC, 
             betas = betas,
             D = D,
             b = b)

res <- gibbs_sim(data)
unlist(res$betas)

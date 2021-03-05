sim_mstate <- function(seed, nsubj = 500) {
  require(splines)
  require(JMbayes)
  require(mstate)
  require(MASS)
  require(Matrix)
  
  set.seed(seed)
  
  N <- nsubj
  
  n <- 20
  
  id <- rep(1:N, each = n)
  
  min.t <- 0.01
  max.t <- 12
  
  time <- replicate(N, c(0, sort(runif(n - 1, min = min.t, max = max.t))), simplify = FALSE)
  time <- do.call(c, time)
  
  Xcov.s <- rnorm(N, mean = 4.763, sd = 2.8)
  Xcov <- rep(Xcov.s, each = n)
  
  Xcov2.s <- sample(c(0, 1), N, replace = T)
  Xcov2 <- rep(Xcov.s, each = n)
  
  simDF <- data.frame("id" = id, "time" = time, "X" = Xcov, "X2" = Xcov2)
  
  X <- model.matrix(~ 1 + time + X, data = simDF)
  Z <- model.matrix(~ 1 + time, data = simDF)
  
  D <- matrix(c(1.35, 4, 4, 2), ncol = 2, nrow = 2)
  D <- as.matrix(Matrix:::nearPD(D)$mat)
  
  b <- MASS:::mvrnorm(n = N, mu = rep(0, ncol(D)), D)
  
  true.betas <- c(-0.482, 0.243, 1.52)
  
  eta.y <- as.vector(X %*% true.betas + rowSums(Z * b[id, ]))
  
  sigma.e <- 1.242
  
  alpha <- c("alpha.01" = 0.8, "alpha.02" = 0.55, "alpha.12" = 1.45, 
             "alpha.sl.01" = 0, "alpha.sl.02" = 0, "alpha.sl.12" = 0, 
             "alpha.ar.01" = 0, "alpha.ar.02" = 0, "alpha.ar.12" = 0)
  
  phi <- c("phi.01" = 12.325, "phi.02" = 12.216, "phi.12" = 10.657)
  
  gammas <- c("(Intercept)1" = -22.25, "X.01" = 1.2,  
              "(Intercept)2" = -18.25, "X.02" = 0.75,
              "(Intercept)3" = -17.25, "X.12" = 0.5)
  
  W <- cbind("(Intercept)1"= rep(1, N), Xcov[seq(1, by = n, N*n)], 
             "(Intercept)2"= rep(1, N), Xcov[seq(1, by = n, N*n)],
             "(Intercept)3"= rep(1, N), Xcov[seq(1, by = n, N*n)])
  
  eta.t1 <- as.vector(W[, c(1, 2), drop = FALSE] %*% gammas[1:2])
  ## transition: 0 -> 2
  eta.t2 <- as.vector(W[, c(3, 4), drop = FALSE] %*% gammas[3:4])
  ## transition: 1 -> 2
  eta.t3 <- as.vector(W[, c(5, 6), drop = FALSE] %*% gammas[5:6])
  
  invS01 <- function(t, u, i) {
    h <- function(s) {
      XX <- cbind(1, s, Xcov[i])
      ZZ <- cbind(1, s)
      XXslope <- cbind(0, rep(1, length(s)))
      ZZslope <- cbind(0, rep(1, length(s)))
      XXarea <- cbind(s, (s^2)/2, s*Xcov[i])
      ZZarea <- cbind(s, (s^2)/2)
      f1 <- as.vector(XX %*% true.betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
      #f1.sl <- as.vector(XXslope %*% true.betas[1:2] + rowSums(ZZslope * b[rep(i, nrow(ZZ)), ]))
      #f1.ar <- as.vector(XXarea %*% true.betas + rowSums(ZZarea * b[rep(i, nrow(ZZ)), ]))
      #exp(log(phi["phi.01"]) + (phi["phi.01"])*log(s) + eta.t1[i] + f1*alpha["alpha.01"] + 
      #      f1.sl*alpha["alpha.sl.01"] + f1.ar*alpha["alpha.ar.01"])
      exp(log(phi["phi.01"]) + (phi["phi.01"])*log(s) + eta.t1[i] + f1*alpha["alpha.01"])
    }
    integrate(h, lower = 0, upper = t, subdivisions = 10000)$value + log(u)
  }
  
  invS02 <- function(t, u, i) {
    h <- function(s) {
      XX <- cbind(1, s, Xcov[i])
      ZZ <- cbind(1, s)
      #XXslope <- cbind(0, rep(1, length(s)))
      #ZZslope <- cbind(0, rep(1, length(s)))
      #XXarea <- cbind(s, (s^2)/2, s*Xcov[i])
      #ZZarea <- cbind(s, (s^2)/2)
      f1 <- as.vector(XX %*% true.betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
      #f1.sl <- as.vector(XXslope %*% true.betas[1:2] + rowSums(ZZslope * b[rep(i, nrow(ZZ)), ]))
      #f1.ar <- as.vector(XXarea %*% true.betas + rowSums(ZZarea * b[rep(i, nrow(ZZ)), ]))
      #exp(log(phi["phi.02"]) + (phi["phi.02"])*log(s) + eta.t2[i] + f1*alpha["alpha.02"] + 
      #      f1.sl*alpha["alpha.sl.02"] + f1.ar*alpha["alpha.ar.02"])
      exp(log(phi["phi.02"]) + (phi["phi.02"])*log(s) + eta.t2[i] + f1*alpha["alpha.02"])
    }
    integrate(h, lower = 0, upper = t, subdivisions = 10000)$value + log(u)
  }
  
  invS12 <- function(t, u, i) {
    h <- function(s) {
      XX <- cbind(1, s, Xcov[i])
      ZZ <- cbind(1, s)
      #XXslope <- cbind(0, rep(1, length(s)))
      #ZZslope <- cbind(0, rep(1, length(s)))
      #XXarea <- cbind(s, (s^2)/2, s*Xcov[i])
      #ZZarea <- cbind(s, (s^2)/2)
      f1 <- as.vector(XX %*% true.betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
      #f1.sl <- as.vector(XXslope %*% true.betas[1:2] + rowSums(ZZslope * b[rep(i, nrow(ZZ)), ]))
      #f1.ar <- as.vector(XXarea %*% true.betas + rowSums(ZZarea * b[rep(i, nrow(ZZ)), ]))
      #exp(log(phi["phi.12"]) + (phi["phi.12"])*log(s) + eta.t3[i] + f1*alpha["alpha.12"] + 
      #      f1.sl*alpha["alpha.sl.12"] + f1.ar*alpha["alpha.ar.12"])
      exp(log(phi["phi.12"]) + (phi["phi.12"])*log(s) + eta.t3[i] + f1*alpha["alpha.12"])
    }
    integrate(h, lower = 0, upper = t, subdivisions = 10000)$value + log(u)
  }
  
  # Probability for each transition
  u01 <- runif(N, 0, 1)
  u02 <- runif(N, 0, 1)
  u12 <- runif(N, 0, 1)
  
  trueT01 <- numeric(N)
  trueT02 <- numeric(N)
  trueT12 <- numeric(N)
  
  mean.Cens <- 9
  C <- runif(N, 0, 2 * mean.Cens)
  
  for (i in 1:N) {
    Root01 <- NULL
    Root02 <- NULL
    Root12 <- NULL
    
    Up <- 50
    tries <- 5
    # Transition 0->1
    Root01 <- try(uniroot(invS01, interval = c(1e-05, Up), u = u01[i], i = i)$root, TRUE)
    while(inherits(Root01, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 200
      Root01 <- try(uniroot(invS01, interval = c(1e-05, Up), u = u01[i], i = i)$root, TRUE)
    }
    trueT01[i] <- if (!inherits(Root01, "try-error")) Root01 else 500
    
    # Transition 0->2
    Root02 <- try(uniroot(invS02, interval = c(1e-05, Up), u = u02[i], i = i)$root, TRUE)
    while(inherits(Root02, "try-error") && tries > 0) {
      tries <- tries - 1
      Up <- Up + 200
      Root02 <- try(uniroot(invS02, interval = c(1e-05, Up), u = u02[i], i = i)$root, TRUE)
    }
    trueT02[i] <- if (!inherits(Root02, "try-error")) Root02 else 500
    
    # Transition 1->2
    if(as.numeric(trueT01[i]) < as.numeric(trueT02[i]) && as.numeric(trueT01[i]) < C[i]) {
      Root12 <- try(uniroot(invS12, interval = c(as.numeric(trueT01[i]), Up), u = u12[i], i = i)$root, TRUE)
      while(inherits(Root12, "try-error") && tries > 0) {
        tries <- tries - 1
        Up <- Up + 200
        Root12 <- try(uniroot(invS12, interval = c(as.numeric(trueT01[i]), Up), u = u12[i], i = i)$root, TRUE)
      }
    } else {Root12 <- 500}
    trueT12[i] <- if (!inherits(Root12, "try-error")) Root12 else 500
    #cat(paste("Subject:", i), "\n")
  }
  
  matsurv <- NULL
  datasurv <- NULL
  for(k in 1:N){
    if (C[k] < min(trueT01[k], trueT02[k])) { # if 0 -> C
      aux1 <- c(k, Xcov.s[k], Xcov2.s[k], C[k], 0, C[k], 0)
      matsurv <- rbind(matsurv, aux1) 
    } else {
      if (trueT02[k] < trueT01[k]) { # if 0 -> 2
        aux1 <- c(k, Xcov.s[k], Xcov2.s[k], trueT02[k], 0, trueT02[k], 1)
        matsurv <- rbind(matsurv, aux1)
      } else {
        if (C[k] < trueT12[k]) { # if 0 -> 1 -> C
          aux1 <- c(k, Xcov.s[k], Xcov2.s[k], trueT01[k], 1, C[k], 0)
          matsurv <- rbind(matsurv, aux1)
        } else { # if 0 -> 1 -> 2
          aux1 <- c(k, Xcov.s[k], Xcov2.s[k], trueT01[k], 1, trueT12[k], 1)
          matsurv <- rbind(matsurv, aux1)
        }
      }
    }
  }
  
  y <- rnorm(N*n, eta.y, sigma.e)
  
  matlongit <- NULL
  aux2 <- NULL
  datalongit <- NULL
  for(k in 1:N){
    n_final <- NULL
    n_final <- if (matsurv[k,4] == 1){
      sum(time[(n*(k-1)+1) : (k*n)] < trueT12[k])
    } else if (matsurv[k,6] == 1){
      sum(time[(n*(k-1)+1) : (k*n)] < trueT02[k])
    } else{
      sum(time[(n*(k-1)+1) : (k*n)] < C[k])          
    }
    aux2 <- matrix(nrow = n_final, ncol = 4,
                   c(rep(k, n_final),
                     y[(n*(k-1)+1) : (n*(k-1) + n_final)],
                     time[(n*(k-1)+1) : (n*(k-1) + n_final)], 
                     Xcov[(n*(k-1)+1) : (n*(k-1) + n_final)]))
    matlongit <- rbind(matlongit, aux2)
  }
  
  
  datasurv <- data.frame(matsurv, row.names = NULL)
  names(datasurv) <- c("id", "X", "X2", "t_State1", "State1", "t_State2", "State2")
  row.names(datasurv) <- as.integer(1:nrow(datasurv))
  
  datalongit <- data.frame(matlongit, row.names = NULL)
  names(datalongit) <- c("id", "Y", "times", "X")
  row.names(datalongit) <- as.integer(1:nrow(datalongit))
  
  simDF$y <- y
  
  tmat <- matrix(NA, 3, 3)
  tmat[1, 2:3] <- 1:2
  tmat[2, 3] <- 3
  dimnames(tmat) <- list(from = c("State0", "State1", "State2"), 
                         to = c("State0", "State1", "State2"))
  
  covs <- c("X", "X2")
  data_mstate <- msprep(time = c(NA, "t_State1", "t_State2"), 
                        status = c(NA, "State1", "State2"), 
                        data = datasurv, 
                        keep = covs,
                        id = "id",  
                        trans = tmat)
  
  data_mstate <- expand.covs(data_mstate, covs, 
                             append = TRUE, longnames = FALSE)
  
  trans1.evnts <- table(data_mstate[data_mstate$trans == 1, "status"])
  trans2.evnts <- table(data_mstate[data_mstate$trans == 2, "status"])
  trans3.evnts <- table(data_mstate[data_mstate$trans == 3, "status"])
  
  range(datalongit$times)
  median(tapply(datalongit$times, datalongit$id, length))
  
  out <- list("params" = list(N = N, n = n, knots = knots, D = D, alpha = alpha, phi = phi, 
                              sigma = sigma.e, gammas = gammas, eta.y = eta.y), 
              "event_summary" = list("trans01" = trans1.evnts, "trans02" = trans3.evnts, 
                                     "trans12" = trans3.evnts), 
              "datasets" = list("datasurv" = datasurv, "datalongit" = datalongit, 
                                "data_mstate" = data_mstate))
  out
}

library("survival")
library("nlme")
library("GLMMadaptive")
library("splines")
#library("Formula")
data("pbc2", package = "JM")
data("pbc2.id", package = "JM")
source(file.path(getwd(), "R/jm.R"))
source(file.path(getwd(), "R/help_functions.R"))
source(file.path(getwd(), "Development/jm/R_to_Cpp.R"))
source(file.path(getwd(), "Development/jm/PBC_data.R"))

fm1 <- lme(log(serBilir) ~ year * sex + I(year^2) + age + prothrombin,
           data = pbc2, random = ~ year | id)
fm2 <- lme(serChol ~ ns(year, 3) + sex + age, data = pbc2, random = ~ year | id,
           na.action = na.exclude)
fm3 <- mixed_model(hepatomegaly ~ sex + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())
fm4 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ 1 | id, family = binomial())
Mixed_objects <- list(fm1, fm2, fm3, fm4)

D_lis <- lapply(Mixed_objects, extract_D)
D <- bdiag(D_lis)
invD <- solve(D)

##########################################################################################
##########################################################################################
dmvnorm <- function (x, mu, Sigma = NULL, invSigma = NULL, log = TRUE,
                     prop = TRUE) {
  if (!is.matrix(x))
    x <- rbind(x)
  if (is.matrix(mu)) {
    if (nrow(mu) != nrow(x))
      stop("incorrect dimensions for 'mu'.")
    p <- ncol(mu)
  }
  else {
    p <- length(mu)
    mu <- rep(mu, each = nrow(x))
  }
  if (is.null(Sigma) && is.null(invSigma))
    stop("'Sigma' or 'invSigma' must be given.")
  if (!is.null(Sigma)) {
    if (is.list(Sigma)) {
      ev <- Sigma$values
      evec <- Sigma$vectors
    } else {
      ed <- eigen(Sigma, symmetric = TRUE)
      ev <- ed$values
      evec <- ed$vectors
    }
    invSigma <- evec %*% (t(evec)/ev)
    if (!prop)
      logdetSigma <- sum(log(ev))
  } else {
    if (!prop)
      logdetSigma <- -determinant(as.matrix(invSigma))$modulus
  }
  ss <- x - mu
  quad <- 0.5 * rowSums((ss %*% invSigma) * ss)
  if (!prop)
    fact <- -0.5 * (p * log(2 * pi) + logdetSigma)
  if (log) {
    if (!prop)
      as.vector(fact - quad)
    else as.vector(-quad)
  } else {
    if (!prop)
      as.vector(exp(fact - quad))
    else as.vector(exp(-quad))
  }
}

cor2cov <- function (R, vars) {
  p <- nrow(R)
  if (length(vars) != p) {
    stop("incorrect length of 'vars' argument.")
  }
  sds <- sqrt(vars)
  sds * R * rep(sds, each = p)
}

ddirichlet <- function (x, alpha, log = FALSE) {
  dirichlet1 <- function(x, alpha) {
    logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    s <- sum((alpha - 1) * log(x))
    exp(sum(s) - logD)
  }
  if (!is.matrix(x))
    if (is.data.frame(x))
      x <- as.matrix(x)
    else x <- t(x)
    if (!is.matrix(alpha))
      alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x),
                      byrow = TRUE)
    if (any(dim(x) != dim(alpha)))
      stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
    pd <- vector(length = nrow(x))
    for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i,
                                                           ])
    pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
    pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
    if (log) return(log(pd)) else return(pd)
}

rdirichlet <- function (n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  x / rowSums(x)
  #exp(log(x) - log(rowSums(x)))
}

## Added half-t distribution
## Trying to implement lkj prior as in stan_mvmer (see http://mc-stan.org/rstanarm/reference/priors.html)

dhalf_t <- function(x, scale = 10, df = 1, logscale = FALSE) {
  if (any(scale < 0, df < 0))
    stop("Scale and df must be positive")
  out <- numeric(length(x))
  out[x <= 0] <- 0
  constant.num <- 2*gamma((df + 1)/2)
  constant.den <- gamma(df/2)*sqrt(df*pi*(scale^2))
  kernel <- (1 + (1 / df) * (x/scale)^2)^(-(df + 1) / 2)
  out[x > 0] <- (constant.num / constant.den) * kernel
  if (logscale == TRUE) {
    return(log(out))
  } else {
    return(out)
  }
}

##########################################################################################
##########################################################################################


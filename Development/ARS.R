library("JMbayes2")
lmeFit <- lme(log(serBilir) ~ poly(year, 3) * sex, data = pbc2,
              random = list(id = pdDiag(form = ~ poly(year, 3))))
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
jointFit <- jm(CoxFit, lmeFit, time_var = "year")

################################################################################
################################################################################
################################################################################

obj <- jointFit
idL_lp <- obj$model_data$idL_lp[[1]]
y <- obj$model_data$y[[1]]
X <- obj$model_data$X[[1]]
Z <- obj$model_data$Z[[1]]
betas <- obj$statistics$Mean$betas1
b <- obj$statistics$Mean$b
sigma <- obj$statistics$Mean$sigmas
D <- matrix(0, ncol(Z), ncol(Z))
D[lower.tri(D, TRUE)] <- obj$statistics$Mean$D
##
id_GK <- rep(1:312, each = 15)
log_Pwk <- obj$model_data$log_Pwk
Time_right <- obj$model_data$Time_right
delta <- obj$model_data$delta
W0_h <- obj$model_data$W0_h
W0_H <- obj$model_data$W0_H
W_h <- obj$model_data$W_h
W_H <- obj$model_data$W_H
X_h <- obj$model_data$X_h[[1]][[1]]
X_H <- obj$model_data$X_H[[1]][[1]]
Z_h <- obj$model_data$Z_h[[1]][[1]]
Z_H <- obj$model_data$Z_H[[1]][[1]]
bs_gammas <- obj$statistics$Mean$bs_gammas
gammas <- obj$statistics$Mean$gammas
alphas <- obj$statistics$Mean$alphas

#####

i <- 33
yi <- y[idL_lp == i]
Xi <- X[idL_lp == i, ]
Zi <- Z[idL_lp == i, ]
deltai <- delta[i]
X_hi <- X_h[i, ]
Z_hi <- Z_h[i, ]
log_Pwki <- log_Pwk[id_GK == i]
W0_Hi <- W0_H[id_GK == i, ]
W_Hi <- W_H[id_GK == i, , drop = FALSE]
X_Hi <- X_H[id_GK == i, ]
Z_Hi <- Z_H[id_GK == i, ]

bi <- b[i, ]

log_post_cond <- function (bi) {
    mui <- c(Xi %*% betas) + c(Zi %*% bi)
    log_pyb <- sum(dnorm(yi, mui, sigma, log = TRUE))
    log_pb <- -0.5 * c(crossprod(bi, solve(D, bi)))
    log_ptb <- deltai * alphas * (c(X_hi %*% betas) + c(Z_hi %*% bi)) -
        sum(exp(log_Pwki + c(W0_Hi %*% bs_gammas + W_Hi %*% gammas) +
                alphas * (c(X_Hi%*% betas) + c(Z_Hi %*% bi))))
    c(log_pyb + log_ptb + log_pb)
}
gradient <- function (bi) {
    c1 <- c(crossprod(Zi, yi - Xi %*% betas - Zi %*% bi)) / sigma^2
    c2 <- alphas * (deltai * Z_hi - colSums(Z_Hi *
        exp(log_Pwki + c(W0_Hi %*% bs_gammas + W_Hi %*% gammas) +
                    alphas * (c(X_Hi%*% betas) + c(Z_Hi %*% bi)))))
    c3 <- - solve(D, bi)
    c1 + c2 + c3
}
gradient2 <- function (bi) JMbayes2:::cd(bi, log_post_cond)

gradient(bi)
gradient2(bi)

hess <- function (bi) JMbayes2:::cd_vec(bi, gradient)

gradient(bi)
hess(bi)

bvals <- seq(-200, 200, len = 500)
yvals <- sapply(bvals, function (bb, pos) {
    bi[pos] <- bb
    log_post_cond(bi)
}, pos = 4)

plot(bvals, yvals, type = "l")



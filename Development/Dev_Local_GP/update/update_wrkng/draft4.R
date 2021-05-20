JMbayes:::update.mvJMbayes
JMbayes:::update.JMbayes

function (object, ...) 
{
  call <- object$call
  if (is.null(call)) 
    stop("need an object with call component.\n")
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras) > 0) {
    nams <- names(extras)
    existing <- !is.na(match(nams, names(call)))
    for (a in names(extras)[existing]) {
      call[[a]] <- extras[[a]]
    }
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
    if (all(nams %in% c("scales", "n.iter", "n.burnin", 
                        "n.adapt", "n.thin"))) {
      call <- c(as.list(call), list(init = extractInits(object)))
      call <- as.call(call)
    }
  }
  else {
    call <- c(as.list(call), list(init = extractInits(object)))
    call <- as.call(call)
  }
  eval(call, parent.frame())
}

i <- 1

x <- joint_model_fit_1
names(x_mcmc$b)
length(x_mcmc[namsX %in% 'b'])
lapply(list, function)


extract_last_iterations <- function(x) {
  x_mcmc <- x$mcmc
  n_chains <- x$control$n_chains
  namsX <- names(x_mcmc)
  last_iter_b <- list(NULL)
  last_iter <- list(NULL)
  for (i in 1:n_chains) {
    last_iter_b[[i]] <- lapply(x_mcmc[namsX %in% 'b'], function(x, i) x[[i]][, , dim(x[[i]])[length(dim(x[[i]]))], drop = FALSE], i = i)
    last_iter[[i]] <- lapply(x_mcmc[!namsX %in% 'b'], function(x, i) x[[i]][nrow(x[[i]]), , drop = FALSE], i = i)
    last_iter[[i]] <- c(last_iter[[i]], last_iter_b[[i]])
  }
  last_iter
}

chk <- extract_last_iteration(x)

chk1 <- last_iter_b[[i]]$b
chk2 <- x$mcmc$b$`1`

chk3 <- last_iter[[i]]$bs_gammas
chk4 <- x$mcmc$bs_gammas$`1`[6000, ]

all.equal(chk3, chk4)
class(chk1)

length(chk)
chk[[2]]$bs_gammas

x$mcmc$bs_gammas$`2`[6000, ]

last_iterations[[1]]$betas <- list(last_iterations[[1]]$betas1)
last_iterations[[1]]$b <- list(last_iterations[[1]]$b)
last_iterations[[1]]$D <- 
initial_values$D
last_iterations[[1]]$
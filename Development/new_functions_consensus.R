# Functions ====================================================================
slicer <- function (n_slices, id_var, data_long, data_surv, seed = 123L) {
  stopifnot(is.numeric(n_slices), length(n_slices) == 1, n_slices >= 1,
            is.character(id_var))
  if (!id_var %in% names(data_long)) {
    stop(paste0("'", id_var, "' not found in data_long."))
  }
  if (!id_var %in% names(data_surv)) {
    stop(paste0("'", id_var, "' not found in data_surv."))
  }
  ids_unq <- unique(c(as.character(data_long[[id_var]]), 
                      as.character(data_surv[[id_var]])))
  if (!is.null(seed)) set.seed(seed)
  grp <- ((seq_along(ids_unq) - 1) %% n_slices) + 1 # assign each ID a group number 1...n_slices in round-robin order (1,2,...,n_slices,1,2,...)
  ids_slc <- split(sample(ids_unq), grp)
  long <- lapply(ids_slc, function (ids) data_long[data_long[[id_var]] %in% ids, ])
  surv <- lapply(ids_slc, function (ids) data_surv[data_surv[[id_var]] %in% ids, ])
  class(long) <- c("sliced_data", class(long))
  class(surv) <- c("sliced_data", class(surv))
  list(long = long, surv = surv)
}


# nlme::lme() / survival::coxph() do not dispatch on the class of the `data` argument.
# We define our own S3 generics that dispatch on `data` (2nd argument),
# so that lme.sliced_data() / coxph.sliced_data() are used when `data` is `sliced_data`,
# while lme.data.frame() / coxph.data.frame() are used for `data.frames`
lme <- function (fixed, data, ...) UseMethod ("lme", data)
mixed_model <- function (fixed, data, ...) UseMethod ("mixed_model", data)
coxph <- function (formula, data, ...) UseMethod ("coxph", data)
jm <- function (Surv_object, Mixed_objects, ...) UseMethod ("jm", Surv_object)

# We cannot implement coxph.data.frame() as a simple forwarder like
# coxph.data.frame <- function (formula, data, ...) {
#   survival::coxph(formula = formula, data = data, ...)
# }
# because the resulting coxph object will typically store the call with
# `object$call$data` as `data`, where `data` is just the argument name of this 
# method, not the original expression supplied by the user (e.g., `pbc2.id`).
# Later when jm() tries to recover the survival dataset via 
# `eval(Surv_object$call$data, parent.frame())`, if `object$call$data` is `data`, 
# jm() may end up evaluating whatever object is named `data` in the environment.
# For this reason, we reconstruct the call using match.call() and evaluate it in the
# parent frame, so object$call$data retains the original expression (e.g., `pbc2.id`).
coxph.default <- function (formula, data, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- quote(survival::coxph)
  eval(mc, parent.frame())
}

lme.default <- function (fixed, data, random, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- quote(nlme::lme)
  eval(mc, parent.frame())
}

mixed_model.default <- function (fixed, random, data, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- quote(GLMMadaptive::mixed_model)
  eval(mc, parent.frame())
}

jm.default <- function (Surv_object, Mixed_objects, ...) {
  mc <- match.call(expand.dots = TRUE)
  mc[[1L]] <- quote(JMbayes2::jm)
  eval(mc, parent.frame())
}

slapply <- function (X, FUN, ..., parallel = TRUE, ncores = NULL,
                    pkgs = NULL, backend = c("auto", "multicore", "snow")) {
  backend <- match.arg(backend)
  
  if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)
  
  if (!parallel || ncores <= 1) return(lapply(X, FUN, ...))
  
  if (backend == "auto") {
    # multicore: parallel::mclapply() uses OS forking (Unix/macOS/Linux only). 
    #            Fast/low overhead because workers are forked copies of the 
    #            current R process. Downside: forking can be fragile if the 
    #            current R session has loaded "native" compiled code (C/C++/Fortran) 
    #            that uses threads internally (e.g., %*%, chol()). In those cases, 
    #            forked workers may hang/crash or behave unpredictably.
    # snow: PSOCK socket cluster, parallel::makeCluster(), starts fresh R worker 
    #       sessions. Slower startup and more data transfer (packages must be 
    #       loaded on workers), but generally the most robust and cross-platform 
    #       option (Windows/Linux/macOS).
    backend <- if (.Platform$OS.type == "windows") "snow" else "multicore"
  }
  
  if (backend == "multicore") {
    return(parallel::mclapply(X, FUN, ..., mc.cores = ncores))
  }
  
  # backend == "snow"
  cl <- parallel::makeCluster(ncores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)
  
  if (!is.null(pkgs)) { # load packages on workers
    pkgs <- as.character(pkgs)
    load_pkgs <- function (pkgs) {
      for (p in pkgs) {
        suppressPackageStartupMessages(library (p, character.only = TRUE))
      }
      NULL
    }
    environment(load_pkgs) <- baseenv()
    parallel::clusterCall(cl, load_pkgs, pkgs = pkgs)
  }
  
  parallel::parLapply(cl, X, FUN, ...)
}

coxph.sliced_data <- function (formula, data, ...,
                               parallel_out = TRUE, cores = NULL) {
  dots <- list(...)
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
  
  fit  <- function (dt, formula, dots) {
    args <- c(list(formula = formula, data = dt), dots)
    do.call(survival::coxph, args)
  }
  
  environment(fit) <- new.env(parent = baseenv()) # Detach fit from the parent call frame so it doesn't accidentally serialize large objects (i.e., the full `data` list) when sent to "snow" workers. print(ls(environment(fit)))
  
  fits <- slapply(X = data, FUN = fit, 
                  formula = formula, dots = dots,
                  parallel = parallel_out, ncores = cores,
                  pkgs = "survival")
  
  class(fits) <- c("sliced_coxph", class(fits))
  fits
}

lme.sliced_data <- function (fixed, data, random, ..., 
                             parallel_out = TRUE, cores = NULL) {
  dots <- list(...)
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
  
  fit <- function (dt, fixed, random, dots) {
    args <- c(list(fixed = fixed, data = dt, random = random), dots)
    do.call(nlme::lme, args)
  }
  
  environment(fit) <- new.env(parent = baseenv()) # Detach fit from the parent call frame so it doesn't accidentally serialize large objects (i.e., the full `data` list) when sent to "snow" workers. print(ls(environment(fit)))
  
  fits <- slapply(X = data, FUN = fit,
                  fixed = fixed, random = random, dots = dots,
                  parallel = parallel_out, ncores = cores,
                  pkgs = "nlme")
  
  class(fits) <- c("sliced_lme", class(fits))
  fits
}

mixed_model.sliced_data <- function (fixed, data, random, ..., 
                                     parallel_out = TRUE, cores = NULL) {
  dots <- list(...)
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
  
  bool1 <- parallel_out && cores > 1L 
  bool2 <- identical(dots$optimizer, "optimParallel") || 
    (is.list(dots$control) && identical(dots$control$optimizer, "optimParallel")) # identical(NULL, "optimParallel") -> FALSE
  if (bool1 && bool2) {
    warning(
      "'optimParallel' requested, but mixed_model.sliced_data() is already parallel across subsamples; ",
      "switching optimizer to 'optim' to avoid nested parallelism. ",
      "To use 'optimParallel', run mixed_model.sliced_data(parallel_out = FALSE).",
      call. = FALSE
    )
    if (identical(dots$optimizer, "optimParallel")) dots$optimizer <- "optim"
    if (is.list(dots$control) && identical(dots$control$optimizer, "optimParallel")) {
      dots$control$optimizer <- "optim"
    }
  }
  
  fit <- function (dt, fixed, random, dots) {
    args <- c(list(fixed = fixed, random = random, data = dt), dots)
    do.call(GLMMadaptive::mixed_model, args)
  }
  
  environment(fit) <- new.env(parent = baseenv()) # Detach fit from the parent call frame so it doesn't accidentally serialize large objects (i.e., the full `data` list) when sent to "snow" workers. print(ls(environment(fit)))
  
  fits <- slapply(X = data, FUN = fit, 
                  fixed = fixed, random = random, dots = dots, 
                  parallel = parallel_out, ncores = cores,
                  pkgs = "GLMMadaptive", backend = "snow")
  
  class(fits) <- c("sliced_MixMod", class(fits))
  fits
}

jm.sliced_coxph <- function (Surv_object, Mixed_objects, time_var, ...,
                             parallel_out = TRUE, cores = NULL) {
  n_slices <- length(Surv_object)
  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)
  
  dots <- list(...)
  ctrl <- dots$control
  if (is.null(ctrl)) ctrl <- list()
  
  # n_chains: control$n_chains > n_chains > default
  if (!is.null(ctrl$n_chains)) {
    n_chains <- ctrl$n_chains
  } else if (!is.null(dots$n_chains)) {
    n_chains <- dots$n_chains
  } else {
    n_chains <- 3L
  }
  
  # inner jm() backend (chains)
  # parallel: control$parallel > parallel > auto
  if (is.null(ctrl$parallel) && !is.null(dots$parallel)) ctrl$parallel <- dots$parallel
  if (is.null(ctrl$parallel)) {
    ctrl$parallel <- if (.Platform$OS.type == "windows") "snow" else "multicore"
  }
  if (identical(ctrl$parallel, "multicore") && .Platform$OS.type == "windows") {
    warning("control$parallel='multicore' is not available on Windows; using 'snow'.",
            call. = FALSE)
    ctrl$parallel <- "snow"
  }
  
  inner_ncores <- max(1L, min(n_chains, cores))
  ctrl$cores <- inner_ncores
  ctrl$n_chains <- n_chains
  
  dots[c("n_chains", "cores", "parallel")] <- NULL # Avoid passing control-like args twice. 
  dots$control <- ctrl
  
  # outer cores (across slices)
  outer_ncores <- min(n_slices, max(1L, floor(cores / inner_ncores)))
  
  # Optional: on Windows, avoid nested PSOCK clusters by default
  # if (.Platform$OS.type == "windows" && inner_ncores > 1L && identical(ctrl$parallel, "snow")) {
  #   outer_ncores <- 1L
  # }
  
  # Detect whether Mixed_objects is:
  # i)  a sliced object: list of lme/MixMod fits (one per slice)
  # ii) a list of sliced objects: list(outcome1_sliced, outcome2_sliced, ...) -> transpose
  if (!inherits(Mixed_objects[[1]], c("lme", "MixMod"))) {
    Mixed_objects <- lapply(seq_len(n_slices), 
                            function (i) lapply(Mixed_objects, `[[`, i))
  }
  
  jobs <- Map(function(S, M) list(Surv_object = S, Mixed_objects = M), # Each worker receives only slice-specific inputs.
               Surv_object, Mixed_objects)
  
  fit <- function (job, time_var, dots) {
    args <- c(list(Surv_object = job$Surv_object,
                   Mixed_objects = job$Mixed_objects,
                   time_var = time_var),
              dots)
    do.call(JMbayes2::jm, args)
  }
  
  environment(fit) <- new.env(parent = baseenv())
  
  fits <- slapply(X = jobs, FUN = fit,
                  time_var = time_var, dots = dots,
                  parallel = parallel_out, ncores = outer_ncores, 
                  pkgs = "JMbayes2", backend = "snow")
  
  class(fits) <- c("sliced_jm", class(fits))
  fits
}

consensus <- function (object, parm,
                      method = c("union", "equal_weight", "var_weight"), 
                      seed = 123L) {
  stopifnot(inherits(object, "sliced_jm"), is.character(parm), length(parm) >= 1)
  method <- match.arg(method)
  get_cons <- function (fits, parm, method) {
    slice_mats <- lapply(fits, function (x) do.call(rbind, x$mcmc[[parm]])) # list of [iter, p] matrices (one per slice)
    if (method == "union") {
      return (list(draws = do.call(rbind, slice_mats), # [iter * slice, p]
                  weights = NULL))
    }
    draws <- simplify2array(slice_mats) # [iter, p, slice]
    draws <- apply(draws, c(2, 3), sample)
    ns <- dim(draws)[3]
    snames <- paste0("S", seq_len(ns)) # slice names
    pnames <- dimnames(draws)[[2]]
    if (method == "equal_weight") {
      w <- matrix(1 / ns, nrow = length(pnames), ncol = ns,
                  dimnames = list(pnames, snames))
      return (list(draws = apply(draws, c(1, 2), mean), # [iter, p]
                  weights = w)) 
    }
    # var_weight
    vars <- apply(draws, c(2, 3), var) # [p, slice]
    w <- 1 / pmax(vars, .Machine$double.eps) # [p, slice]
    w <- sweep(w, 1, rowSums(w), "/") # normalize 
    dimnames(w) <- list(pnames, snames)
    w_draws <- sweep(draws, c(2, 3), w, "*") # weighted[iter, p, slice] <- draws[iter, p, slice] * w[p, slice]
    list(draws = apply(w_draws, c(1, 2), sum), weights = w)
  }
  if (!is.null(seed)) set.seed(seed)
  cons_out <- lapply(parm, function (p) get_cons(object, p, method))
  names(cons_out) <- parm
  cons_draws <- lapply(cons_out, `[[`, "draws")
  cons_wts   <- lapply(cons_out, `[[`, "weights")
  n_draws <- vapply(cons_draws, nrow, integer(1))
  summarise_draws <- function (draws_mat) {
    cbind(
      Mean   = colMeans(draws_mat),
      StDev  = apply(draws_mat, 2, sd),
     `2.5%`  = apply(draws_mat, 2, quantile, probs = 0.025, names = FALSE),
     `97.5%` = apply(draws_mat, 2, quantile, probs = 0.975, names = FALSE),
      P      = apply(draws_mat, 2, function (x) 2 * min(mean(x > 0), mean(x < 0)))
      )
  }
  cons_sum <- lapply(cons_draws, summarise_draws)
  res <- list(method = method, parm = parm, n_splits = length(object), 
              draws = cons_draws, weights = cons_wts, n_draws = n_draws,
              summary = cons_sum, seed = seed)
  class(res) <- "consensus_jm"
  res
}

print.consensus_jm <- function (x, digits = 4) {
  cat("\nConsensus summary (", x$n_splits, " splits)\n", sep = "")
  if (identical(x$method, "union")) {
    cat("Method (union): concatenated draws across splits (no averaging).\n")
  } else if (identical(x$method, "equal_weight")) {
    cat("Method (equal_weight): iteration-wise simple average across splits.\n")
  } else if (identical(x$method, "var_weight")) {
    cat("Method (var_weight): iteration-wise weighted average across splits using inverse-variance weights.\n")
  }
  for (nm in names(x$summary)) {
    cat("\nParameter block: ", nm, 
        " (", x$n_draws[[nm]], " draws)\n", sep = "")
    tab <- x$summary[[nm]]
    w <- x$weights[[nm]]
    if (!is.null(w)) {
      w <- w[rownames(tab), , drop = FALSE]
      colnames(w) <- paste0("w_", colnames(w))
      tab <- cbind(tab, w)
    }
    print(round(tab, digits))
  }
  invisible(x)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Examples =====================================================================
## Default joint model =========================================================
library(JMbayes2)
tstart1 <- Sys.time()
{
  lme1 <- lme(log(serBilir) ~ year * sex, data = pbc2, random = ~ year | id)
  mxm1 <- mixed_model(ascites ~ year + sex, data = pbc2,
                      random = ~ year | id, family = binomial())
  cox1 <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
  jm1 <- jm(cox1, list(lme1, mxm1), time_var = "year")
}
tstop1 <- Sys.time()
difftime(tstop1, tstart1)
## Sliced joint model ==========================================================
### Example 1 ==================================================================
tstart2 <- Sys.time()
{
  data <- slicer(n_slices = 2, id_var = "id", data_long = pbc2, 
                 data_surv = pbc2.id, seed = 123L)
  lme2 <- lme(log(serBilir) ~ year * sex, data = data$long, random = ~ year | id)
  mxm2 <- mixed_model(ascites ~ year + sex, data = data$long,
                      random = ~ year | id, family = binomial(),
                      control = list(optimizer = "optimParallel"))
  cox2 <- coxph(Surv(years, status2) ~ sex, data = data$surv)
  jm2 <- jm(cox2, list(lme2, mxm2), time_var = "year")
}
tstop2 <- Sys.time()
difftime(tstop2, tstart2)
cns_vw <- consensus(jm2, parm = c("gammas", "alphas", "betas1", "betas2"), 
                    method = "var_weight")
cns_eq <- consensus(jm2, parm = c("gammas", "alphas", "betas1", "betas2"), 
                    method = "equal_weight")
cns_un <- consensus(jm2, parm = c("gammas", "alphas", "betas1", "betas2"), 
                    method = "union")
cns_vw
cns_vw$summary$gammas # Mean, StDev, 2.5% , 97.5%, P
cns_vw$draws$gammas # consensus draws
### Example 2 ==================================================================


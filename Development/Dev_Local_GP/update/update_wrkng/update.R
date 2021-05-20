

fm1 <- lme(fixed = log(serBilir) ~ year * sex + I(year^2) +
             age + prothrombin, random =  ~ year | id, data = pbc2)

# [2] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [3] The basic joint model is fitted using a call to jm() i.e.,
joint_model_fit_1 <- jm(fCox1, fm1, time_var = "year", n_burnin = 525)
object <- joint_model_fit_1



update <- function(object, ...) {
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
    if (nams %in% c("n_iter")) {
      call[['n_burnin']] <- 0
      call <- c(as.list(call), list(last_iterations = extract_last_iterations(object)))
      call <- as.call(call)
    }
  } else {
    call <- as.call(call)
  }
  eval(call, parent.frame())
}

chk <- update(object, n_iter = 1)
i <- 1

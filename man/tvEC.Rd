\name{tvEC}
\alias{tvEC}
\alias{tvEC.default}
\alias{tvEC.jm}
\alias{tvEC.coxph}

\title{
Expected Cost for Dynamic Predictions under Left-Truncated and Right-Censored Data
}
\description{
Computes estimates of the expected cost for a given decision threshold (\code{decision_threshold}) at a specified horizon time (\code{Thoriz}), using longitudinal information available up to a given starting time point (\code{Tstart}). The function supports joint models, Cox proportional hazards models and predictions obtained from any method. It provides inverse probability weighting and model-based estimation methods.
}
\usage{
tvEC(object, \dots)

\method{tvEC}{default}(object, Time, event, LT_Time = NULL, Tstart,
     decision_threshold = NULL, Thoriz = NULL, Dt = NULL,
     type_weights = c("model-based", "IPCW"),
     preds_MB = NULL, name_model = NULL, \dots)

\method{tvEC}{jm}(object, newdata, Tstart,
     decision_threshold = NULL, Thoriz = NULL, Dt = NULL,
     type_weights = c("model-based", "IPCW"), name_model = NULL, \dots)

\method{tvEC}{coxph}(object, newdata, Tstart,
     decision_threshold = NULL, Thoriz = NULL, Dt = NULL,
     type_weights = c("model-based", "IPCW"), name_model = NULL, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jm} or \code{coxph}, or a numeric vector containing the dynamic predictions to be evaluated.}

  \item{newdata}{a data.frame that contains the longitudinal and covariate information for the subjects for which prediction of survival probabilities is required. If a \code{jm} object is used, the names of the variables in this data.frame must be the same as in the data.frames that were used to fit the linear mixed effects and the event process model that were supplied as the two first argument of \code{\link{jm}}. If a \code{coxph} is used it only applies to the event process. If the default method is used, the argument is not needed.}

  \item{Tstart}{numeric scalar denoting the time point up to which longitudinal information is to be used to derive predictions.}

  \item{decision_threshold}{numeric scalar between 0 and 1 denoting the decision threshold used in the expected cost.}

  \item{Time, event, LT_Time}{numeric vectors that must be specified when using the default method to compute the metric by directly using a vector of predictions, observed times (\code{Time}), the event indicator (\code{event}) and (possibly) left-truncated times (\code{LT_Time}).}

  \item{Thoriz}{numeric scalar denoting the time point for which a prediction of the survival status is of interest; \code{Thoriz} must be later than \code{Tstart} and either \code{Dt} or \code{Thoriz} must be specified. If \code{Thoriz} is \code{NULL} is set equal to \code{Tstart + Dt}.}

  \item{Dt}{numeric scalar denoting the length of the time interval of prediction; either \code{Dt} or \code{Thoriz} must be specified.}

  \item{type_weights}{character string denoting the type of weights to use to account for censorting. Options are model-based (default) and inverse probability weighting (using the Kaplan-Meier estimate of the censoring distribution, and adjusting for left truncation if necessary).}

  \item{preds_MB}{a numeric vector containing the model-based weights to be used for the individuals censored within (\code{Tstart}, \code{Thoriz}]. It must be specified when using the default method and the \code{model-based} estimator. The vector must have length equal to the number of individuals considered and must include non-missing elements only for the censored individuals.}

  \item{name_model}{character string denoting the name of the model or object used for prediction evaluation; default is \code{NULL}, meaning that the name of the object in the environment will be used.}

  \item{\dots}{additional arguments passed to \code{predict.jm()}.}
}
\value{
A list of class \code{tvMetric} with components:
  \item{metric}{a numeric scalar denoting the estimated expected cost.}
  \item{metric_type}{the character string \code{'EC'}.}
  \item{nr}{a numeric scalar denoting the number of subjects at risk at time \code{Tstart}.}
  \item{nint}{a numeric scalar denoting the number of subjects with an event in (\code{Tstart}, \code{Thoriz}].}
  \item{ncens}{a numeric scalar denoting the number of subjects with a censored time in (\code{Tstart}, \code{Thoriz}].}
  \item{Tstart}{a copy of the \code{Tstart} argument.}
  \item{Thoriz}{a copy of the \code{Thoriz} argument.}
  \item{type_weights}{a copy of the \code{type_weights} argument.}
  \item{metric_info}{a copy of the \code{decision_threshold} argument.}
  \item{LT_MS}{a logical that is \code{TRUE} if left-truncated times have been used.}
  \item{classObject}{the class of \code{object}.}
  \item{nameObject}{the name of the model; either a copy of the \code{name_model} argument, or the name of \code{object}.}
}
\references{
Van Claster, B., et al. (2025). Evaluation of performance measures in predictive artificial
intelligence models to support medical decisions: overview and guidance. \emph{The Lancet Digital Health}.
}
\author{
Arnau Garcia-Fernandez \email{a.garciafernandez@erasmusmc.nl}, Pedro Miranda-Afonso \email{p.mirandaafonso@erasmusmc.nl} and Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{decision_curves}}, \code{\link{tvIWEC}}, \code{\link{tvNB}}, \code{\link{predict}}, \code{\link{jm}}
}
\examples{
\donttest{
# We fit a multivariate joint model
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
fm1 <- lme(log(serBilir) ~ ns(year, 3) * sex, data = pbc2,
           random = ~ ns(year, 3) | id, control = lmeControl(opt = 'optim'))
fm2 <- lme(prothrombin ~ ns(year, 2) * sex, data = pbc2,
           random = ~ ns(year, 2) | id, control = lmeControl(opt = 'optim'))

jointFit <- jm(CoxFit, list(fm1, fm2), time_var = "year", n_chains = 1L)

# Expected cost with decision threshold 0.2, prediction horizon 8
# and landmark time 5
tvEC(jointFit, newdata = pbc2, Tstart = 5, Dt = 3, decision_threshold = 0.2,
     cores = 1L)

# Considering the Cox model and IPW approach for estimation:
tvEC(CoxFit, newdata = pbc2.id, Tstart = 5, Dt = 3, decision_threshold = 0.2,
     type_weights = "IPCW")

# Considering the vectors with predictions, observed times, event indicator and
# left-truncated times. We consider the time window (5,8]. We use the default method:
pbc2$lt <- unlist(by(pbc2, pbc2$id, function(x) rep(runif(1, 0, x$years), nrow(x))))
pbc2.id$lt <- runif(nrow(pbc2.id), 0, pbc2.id$years)
sfit0 <- summary(survfit(CoxFit, newdata = pbc2.id), times = 5)
sfit <- summary(survfit(CoxFit, newdata = pbc2.id), times = 8)
predictions <- 1.0 - (as.matrix(sfit$surv)[1L, ] / as.matrix(sfit0$surv)[1L, ])
tvEC(predictions, pbc2.id$years, pbc2.id$event, pbc2.id$lt, Tstart = 5, Thoriz = 8,
     type_weights = "IPCW", decision_threshold = 0.25)

# Now, we consider the default method, but using the model-based estimator:
ind_censoring <- pbc2.id$years > 5 & pbc2.id$lt < 5 & pbc2.id$years < 8 & pbc2.id$event == 0
sfit2_8 <- summary(survfit(CoxFit, newdata = pbc2.id[ind_censoring, ]),
                 times = 8)
sfit_cens <- summary(survfit(CoxFit, newdata = pbc2.id[ind_censoring, ]),
                     times = pbc2.id$years[ind_censoring])
preds_MB <- rep(NA, length(pbc2.id$years))
preds_MB[ind_censoring] <- 1.0 - (as.matrix(sfit2_8$surv)[1L, ] / as.matrix(sfit_cens$surv)[1L, ])
tvEC(predictions, Time = pbc2.id$years, event = pbc2.id$event, LT_Time = pbc2.id$lt,
     Tstart = 5, Dt = 3, type_weights = "model-based", preds_MB = preds_MB,
     decision_threshold = 0.3)
}
}


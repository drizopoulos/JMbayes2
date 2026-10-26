\name{tvBrier}
\alias{tvBrier}
\alias{tvBrier.jm}
\alias{tvBrier.coxph}

\title{
Brier Score for Dynamic Predictions under Left-Truncated and Right-Censored Data
}
\description{
Computes estimates of the Brier score at a specified horizon time (\code{Thoriz}), using longitudinal information available up to a given starting time point (\code{Tstart}). The function supports both joint models and Cox proportional hazards models, and provides inverse probability weighting and model-based estimation methods.
}
\usage{
tvBrier(object, newdata, Tstart, \dots)

\method{tvBrier}{jm}(object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
          integrated = FALSE, type_weights = c("model-based", "IPCW"),
          model_weights = NULL, eventData_fun = NULL,
          parallel = c("snow", "multicore"),
          cores = parallelly::availableCores(omit = 1L), ...)
      
\method{tvBrier}{coxph}(object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
          integrated = FALSE,
          type_weights = c("model-based", "IPCW"), ...)
}

\arguments{
  \item{object}{an object inheriting from class \code{jm} or \code{coxph}. It can also be a library of joint models.}

  \item{newdata}{a data.frame that contains the longitudinal and covariate information for the subjects for which prediction of survival probabilities is required. If a \code{jm} object is used, the names of the variables in this data.frame must be the same as in the data.frames that were used to fit the linear mixed effects and the event process model that were supplied as the two first argument of \code{\link{jm}}. If a \code{coxph} is used it only applies to the event process.}
  
  \item{Tstart}{numeric scalar denoting the time point up to which longitudinal information is to be used to derive predictions.}
  
  \item{Thoriz}{numeric scalar denoting the time point for which a prediction of the survival status is of interest; \code{Thoriz} must be later than \code{Tstart} and either \code{Dt} or \code{Thoriz} must be specified. If \code{Thoriz} is \code{NULL} is set equal to \code{Tstart + Dt}.}
  
  \item{Dt}{numeric scalar denoting the length of the time interval of prediction; either \code{Dt} or \code{Thoriz} must be specified.}
  
  \item{integrated}{logical; if \code{TRUE} the integrated Brier score is calculated. Integral is computed with respect the time (see below for details).}
  
  \item{type_weights}{character string denoting the type of weights to use to account for censorting. Options are model-based (default) and inverse probability weighting (using the Kaplan-Meier estimate of the censoring distribution, and adjusting for left truncation if necessary).}
  
  \item{model_weights}{a numeric vector of weights to combine predictions when \code{object} is a list of joint models of class \code{"jmList"}.}
  
  \item{eventData_fun}{a function that takes as input the \code{newdata} and produces the dataset used for the event process model. This is useful when, for example, the event process model contains other time-varying covariates.
        It is important that this function does not alter the ordering of the subjects in \code{newdata}.}

  \item{parallel}{character string; what type of parallel computing to use.}
  
  \item{cores}{integer denoting the number of cores to be used when a library of joint models has been provided in
      \code{object}. If \code{cores = 1}, no parallel computing is used.}
      
  \item{\dots}{additional arguments passed to \code{predict.jm()}.}
}
\details{
The integrated Brier score allows us to summarize the predictive accuracy in the interval of interest \eqn{(t,t+s]} by computing the integral \eqn{iBS(t+s,t)=\frac{1}{s}\int^{t+s}_t BS(z,t) dz}. The integral is approximated numerically via Simpson's rule.
}
\value{
A list of class \code{tvBrier} with components:
  \item{Brier}{a numeric scalar denoting the estimated Brier score.}
  \item{nr}{a numeric scalar denoting the number of subjects at risk at time \code{Tstart}.}
  \item{nint}{a numeric scalar denoting the number of subjects with an event in (\code{Tstart}, \code{Thoriz}].}
  \item{ncens}{a numeric scalar denoting the number of subjects with a censored time in (\code{Tstart}, \code{Thoriz}].}
  \item{Tstart}{a copy of the \code{Tstart} argument.}
  \item{Thoriz}{a copy of the \code{Thoriz} argument.}
  \item{integrated}{a copy of the \code{integrated} argument.}
  \item{type_weights}{a copy of the \code{type_weights} argument.}
  \item{LT_MS}{a logical that is \code{TRUE} if left-truncated times have been used.}
  \item{classObject}{the class of \code{object}.}
  \item{nameObject}{the name of \code{object}.}
}
\references{
Rizopoulos, D. (2016). The R package JMbayes for fitting joint models for longitudinal and
time-to-event data using MCMC. \emph{Journal of Statistical Software} \bold{72(7)}, 1--45.
doi:10.18637/jss.v072.i07.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in joint models for
longitudinal and time-to-event data. \emph{Biometrics} \bold{67}, 819--829.

Rizopoulos, D., Molenberghs, G. and Lesaffre, E.M.E.H. (2017). Dynamic predictions with 
time-dependent covariates in survival analysis using joint modeling and landmarking. 
\emph{Biometrical Journal} \bold{59}, 1261--1276.
}
\author{
Arnau Garcia-Fernandez \email{a.garciafernandez@erasmusmc.nl}, Pedro Miranda-Afonso \email{p.mirandaafonso@erasmusmc.nl} and Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}


\seealso{
\code{\link{predict}}, \code{\link{jm}}
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

# We use landmark time 3 and prediction horizon 4.
tvBrier(jointFit, newdata = pbc2, Tstart = 3, Dt = 4)

# Considering the Cox model and the IPW approach for the estimation
tvBrier(CoxFit, newdata = pbc2.id, Tstart = 6, Dt = 4, type_weights = "IPCW")
}
}
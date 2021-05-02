\name{Accuracy Measures}

\alias{tvAUC}
\alias{tvROC}
\alias{tvAUC.jm}
\alias{tvROC.jm}


\title{Time-Dependent ROCs and AUCs for Joint Models}

\description{
  Using the available longitudinal information up to a starting time point, this function computes an estimate of the ROC and the AUC at a horizon time point based on joint models.
}

\usage{
tvAUC(object, newdata, Tstart, \dots)

\method{tvAUC}{jm}(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, \dots)

tvROC(object, newdata, Tstart, \dots)

\method{tvROC}{jm}(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, \dots)
}

\arguments{
\item{object}{an object inheriting from class \code{jm}.}

\item{newdata}{a data.frame that contains the longitudinal and covariate information for the subjects for which prediction of survival probabilities is required. The names of the variables in this data.frame must be the same as in the data.frames that were used to fit the linear mixed effects and the event process model that were supplied as the two first argument of \code{\link{jm}}.}

\item{Tstart}{numeric scalar denoting the time point up to which longitudinal information is to be used to derive predictions.}

\item{Thoriz}{numeric scalar denoting the time point for which a prediction of the survival status is of interest; \code{Thoriz} must be later than \code{Tstart} and either \code{Dt} or \code{Thoriz} must be specified. If \code{Thoriz} is \code{NULL} is set equal to \code{Tstart + Dt}.}

\item{Dt}{numeric scalar denoting the length of the time interval of prediction; either \code{Dt} or \code{Thoriz} must be specified.}

\item{\dots}{additional arguments; currently none is used.}
}

\value{
  A list of class \code{tvAUC} with components:
  \item{auc}{a numeric scalar denoting the estimated prediction error.}
  \item{Tstart}{a copy of the \code{Tstart} argument.}
  \item{Thoriz}{a copy of the \code{Thoriz} argument.}
  \item{nr}{a numeric scalar denoting the number of subjects at risk at time \code{Tstart}.}
  \item{classObject}{the class of \code{object}.}
  \item{nameObject}{the name of \code{object}.}

  A list of class \code{tvROC} with components:
  \item{TP, FP, nTP, nFN, nTN, qSN, qSP, qOverall}{accuracy indexes.}
  \item{F1score, Youden}{numeric scalars with the optimal cut-point using the F1 score and the Youden index.}
  \item{thr}{numeric vector of thresholds.}
  \item{Tstart}{a copy of the \code{Tstart} argument.}
  \item{Thoriz}{a copy of the \code{Thoriz} argument.}
  \item{nr}{a numeric scalar denoting the number of subjects at risk at time \code{Tstart}.}
  \item{classObject}{the class of \code{object}.}
  \item{nameObject}{the name of \code{object}.}
}

\references{
Antolini, L., Boracchi, P., and Biganzoli, E. (2005). A time-dependent discrimination index
for survival data. \emph{Statistics in Medicine} \bold{24}, 3927--3944.

Harrell, F., Kerry, L. and Mark, D. (1996). Multivariable prognostic models: issues in
developing models, evaluating assumptions and adequacy, and measuring and reducing errors.
\emph{Statistics in Medicine} \bold{15}, 361--387.

Heagerty, P. and Zheng, Y. (2005). Survival model predictive accuracy and ROC curves.
\emph{Biometrics} \bold{61}, 92--105.

Rizopoulos, D. (2016). The R package JMbayes for fitting joint models for longitudinal and
time-to-event data using MCMC. \emph{Journal of Statistical Software} \bold{72(7)}, 1--45.
doi:10.18637/jss.v072.i07.

Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Rizopoulos, D. (2011). Dynamic predictions and prospective accuracy in joint models for
longitudinal and time-to-event data. \emph{Biometrics} \bold{67}, 819--829.

Rizopoulos, D., Molenberghs, G. and Lesaffre, E.M.E.H. (2017). Dynamic predictions with time‐dependent covariates in survival analysis using joint modeling and landmarking. \emph{Biometrical Journal} \bold{59}, 1261--1276.
}

\author{Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}}


\seealso{\code{\link{predict}}, \code{\link{jm}}}

\examples{
\donttest{
# We fit a multivariate joint model
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)
fm1 <- lme(log(serBilir) ~ ns(year, 3) * sex, data = pbc2,
           random = ~ ns(year, 3) | id, control = lmeControl(opt = 'optim'))
fm2 <- lme(prothrombin ~ ns(year, 2) * sex, data = pbc2,
           random = ~ ns(year, 2) | id, control = lmeControl(opt = 'optim'))
fm3 <- mixed_model(ascites ~ year * sex, data = pbc2,
                   random = ~ year | id, family = binomial())

jointFit <- jm(CoxFit, list(fm1, fm2, fm3), time_var = "year", n_chains = 1L)

roc <- tvROC(jointFit, newdata = pbc2, Tstart = 4, Dt = 3, cores = 1L)
roc
plot(roc)
}
}

\keyword{methods}

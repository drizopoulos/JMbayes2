\name{tvROC}
\alias{tvROC}
\alias{tvROC.jm}
\alias{tvROC.coxph}

\title{
Receiver Operating Characteristic Curve (ROC) for Dynamic Predictions under Left-Truncated and Right-Censored Data
}
\description{
Computes estimates of the ROC curve at a specified horizon time (\code{Thoriz}), using longitudinal information available up to a given starting time point (\code{Tstart}). The function supports both joint models and Cox proportional hazards models, and provides inverse probability weighting and model-based estimation methods.
}
\usage{
tvROC(object, newdata, Tstart, \dots)

\method{tvROC}{jm}(object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
          type_weights = c("model-based", "IPCW"), ...)
      
\method{tvROC}{coxph}(object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
          type_weights = c("model-based", "IPCW"), ...)
}

\arguments{
  \item{object}{an object inheriting from class \code{jm} or \code{coxph}.}

  \item{newdata}{a data.frame that contains the longitudinal and covariate information for the subjects for which prediction of survival probabilities is required. If a \code{jm} object is used, the names of the variables in this data.frame must be the same as in the data.frames that were used to fit the linear mixed effects and the event process model that were supplied as the two first argument of \code{\link{jm}}. If a \code{coxph} is used it only applies to the event process.}
  
  \item{Tstart}{numeric scalar denoting the time point up to which longitudinal information is to be used to derive predictions.}
  
  \item{Thoriz}{numeric scalar denoting the time point for which a prediction of the survival status is of interest; \code{Thoriz} must be later than \code{Tstart} and either \code{Dt} or \code{Thoriz} must be specified. If \code{Thoriz} is \code{NULL} is set equal to \code{Tstart + Dt}.}
  
  \item{Dt}{numeric scalar denoting the length of the time interval of prediction; either \code{Dt} or \code{Thoriz} must be specified.}
  
  \item{type_weights}{character string denoting the type of weights to use to account for censorting. Options are model-based (default) and inverse probability weighting (using the Kaplan-Meier estimate of the censoring distribution, and adjusting for left truncation if necessary).}
      
  \item{\dots}{additional arguments passed to \code{predict.jm()}.}
}
\value{
A list of class \code{tvROC} with components:
  \item{TP, FP, nTP, nFN, nTN, qSN, qSP, qOverall}{accuracy indexes.}
  \item{F1score, Youden}{numeric scalars with the optimal cut-point using the F1 score and the Youden index.}
  \item{thr}{numeric vector of thresholds.}
  \item{Tstart}{a copy of the \code{Tstart} argument.}
  \item{Thoriz}{a copy of the \code{Thoriz} argument.}
  \item{nr}{a numeric scalar denoting the number of subjects at risk at time \code{Tstart}.}
  \item{classObject}{the class of \code{object}.}
  \item{type_weights}{a copy of the \code{type_weights} argument.}
  \item{LT_MS}{a logical that is \code{TRUE} if left-truncated times have been used.}
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

# We use landmark time 4 and prediction horizon 3.
roc <- tvROC(jointFit, newdata = pbc2, Tstart = 4, Dt = 3, cores = 1L)
roc

# We plot the ROC curve
plot(roc, legend = TRUE, optimal_cutoff = "Youden")
}
}
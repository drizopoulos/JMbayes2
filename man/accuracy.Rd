\name{Accuracy Measures}

\alias{tvAUC}
\alias{tvROC}
\alias{tvAUC.jm}
\alias{tvAUC.tvROC}
\alias{tvROC.jm}
\alias{calibration_plot}
\alias{calibration_metrics}
\alias{tvBrier}
\alias{tvEPCE}
\alias{create_folds}


\title{Time-Dependent Predictive Accuracy Measures for Joint Models}

\description{
  Using the available longitudinal information up to a starting time point, these functions compute estimates of the ROC curve and the AUC, the Brier score and expected predictive cross-entropy at a horizon time point based on joint models.
}

\usage{
tvROC(object, newdata, Tstart, \dots)

\method{tvROC}{jm}(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, type_weights = c("model-based", "IPCW"), \dots)

tvAUC(object, newdata, Tstart, \dots)

\method{tvAUC}{jm}(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, type_weights = c("model-based", "IPCW"), \dots)

\method{tvAUC}{tvROC}(object, \dots)

calibration_plot(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, df_ns = 3, plot = TRUE,
    col = "red", lty = 1, lwd = 1,
    add_CI = TRUE, col_CI = "lightgrey",
    add_density = TRUE, col_dens = "grey",
    xlab = "Predicted Probabilities",
    ylab = "Observed Probabilities", main = "", \dots)

calibration_metrics(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, df_ns = 3, \dots)

tvBrier(object, newdata, Tstart, Thoriz = NULL, Dt = NULL,
    integrated = FALSE, type_weights = c("model-based", "IPCW"),
    model_weights = NULL, eventData_fun = NULL,
    parallel = c("snow", "multicore"),
    cores = parallelly::availableCores(omit = 1L), \dots)

tvEPCE(object, newdata, Tstart, Thoriz = NULL, Dt = NULL, eps = 0.001,
    model_weights = NULL, eventData_fun = NULL,
    parallel = c("snow", "multicore"),
    cores = parallelly::availableCores(omit = 1L), \dots)

create_folds(data, V = 5, id_var = "id", strata = NULL, seed = 123L)
}

\arguments{
\item{object}{an object inheriting from class \code{jm}, except for \code{tvAUC.tvROC()} where this is an object of class \code{tvROC}. For \code{tvBrier()} and \code{tvEPCE()} it can also be a library of joint models.}

\item{newdata}{a data.frame that contains the longitudinal and covariate information for the subjects for which prediction of survival probabilities is required. The names of the variables in this data.frame must be the same as in the data.frames that were used to fit the linear mixed effects and the event process model that were supplied as the two first argument of \code{\link{jm}}.}

\item{Tstart}{numeric scalar denoting the time point up to which longitudinal information is to be used to derive predictions.}

\item{Thoriz}{numeric scalar denoting the time point for which a prediction of the survival status is of interest; \code{Thoriz} must be later than \code{Tstart} and either \code{Dt} or \code{Thoriz} must be specified. If \code{Thoriz} is \code{NULL} is set equal to \code{Tstart + Dt}.}

\item{Dt}{numeric scalar denoting the length of the time interval of prediction; either \code{Dt} or \code{Thoriz} must be specified.}

\item{integrated}{logical; if \code{TRUE} the integrated Brier score is calculated.}

\item{type_weights}{character string denoting the type of weights to use to account for censorting. Options are model-based (default) and inverse probability of censoring weighting (using the Kaplan-Meier estimate of the censoring distribution).}

\item{eps}{numeric scalar used in the approximation of the hazard function.}

\item{model_weights}{a numeric vector of weights to combine predictions when \code{object} is a list of joint models of class \code{"jmList"}.}

\item{eventData_fun}{a function that takes as input the \code{newdata} and produces the dataset used for the event process model. This is useful when, for example, the event process model contains other time-varying covariates.
It is important that this function does not alter the ordering of the subjects in \code{newdata}.}

\item{parallel}{character string; what type of parallel computing to use.}

\item{cores}{integer denoting the number of cores to be used when a library of joint models has been provided in
    \code{object}. If \code{cores = 1}, no parallel computing is used.}

\item{df_ns}{the degrees of freedom for the natural cubic spline of the cloglog transformation of the
predicted probabilities used in the Cox model that assess calibration.}

\item{plot}{logical; should a plot be produced. If \code{FALSE}, a list is returned with the observed and predicted probabilities.}

\item{add_density}{logical; should the kernal density estimation of the predicted probabilities be superimposed in the calibration plot.}

\item{col, lwd, lty, col_dens, xlab, ylab, main}{graphical parameters.}

\item{data}{the data.frame to split in folds.}

\item{V}{numeric scalar denoting the number of folds.}

\item{id_var}{character string denoting the name of the subject id variable in \code{data}.}

\item{strata}{character vector with the names of stratifying variables.}

\item{seed}{integer denoting the seed.}

\item{\dots}{additional arguments passed to \code{predict.jm()}.}
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

Commenges, D., Liquet, B., and Proust-Lima, C. (2012). Choice of prognostic estimators
in joint models by estimating differences of expected conditional Kullback-Leibler risks.
\emph{Biometrics} \bold{68}, 380--387.

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

Rizopoulos, D., Molenberghs, G. and Lesaffre, E.M.E.H. (2017). Dynamic predictions with time-dependent covariates in survival analysis using joint modeling and landmarking. \emph{Biometrical Journal} \bold{59}, 1261--1276.
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
tvAUC(roc)
plot(roc, legend = TRUE, optimal_cutoff = "Youden")
}
}

\keyword{methods}

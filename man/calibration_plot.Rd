\name{calibration_plot}
\alias{calibration_plot}
\alias{calibration_plot.jm}
\alias{calibration_plot.coxph}
\alias{calibration_metrics}

\title{
Calibration Plot for Dynamic Predictions under Right-Censored Data
}
\description{
Displays a calibration plot at a specified horizon time (\code{Thoriz}), using longitudinal information available up to a given starting time point (\code{Tstart}). The function supports both joint models and Cox proportional hazards models. The \code{calibration_metrics()} function calculates metrics for the accuracy of predictions: The ICI is the mean absolute difference between the observed and predicted probabilities, E50 is the median absolute difference, and E90 is the 90\% percentile of the absolute differences.
}
\usage{
calibration_plot(object, newdata, Tstart, \dots)

\method{calibration_plot}{jm}(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, df_ns = NULL, plot = TRUE,
    col = "red", lty = 1, lwd = 1,
    add_CI = TRUE, col_CI = "lightgrey",
    add_density = TRUE, col_dens = "grey",
    xlab = "Predicted Probabilities",
    ylab = "Observed Probabilities", main = "", \dots)
      
\method{calibration_plot}{coxph}(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, df_ns = NULL, plot = TRUE,
    col = "red", lty = 1, lwd = 1,
    add_CI = TRUE, col_CI = "lightgrey",
    add_density = TRUE, col_dens = "grey",
    xlab = "Predicted Probabilities",
    ylab = "Observed Probabilities", main = "", \dots)
          
calibration_metrics(object, newdata, Tstart, Thoriz = NULL,
    Dt = NULL, df_ns = NULL, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{jm} or \code{coxph}.}

  \item{newdata}{a data.frame that contains the longitudinal and covariate information for the subjects for which prediction of survival probabilities is required. If a \code{jm} object is used, the names of the variables in this data.frame must be the same as in the data.frames that were used to fit the linear mixed effects and the event process model that were supplied as the two first argument of \code{\link{jm}}. If a \code{coxph} is used it only applies to the event process.}
  
  \item{Tstart}{numeric scalar denoting the time point up to which longitudinal information is to be used to derive predictions.}
  
  \item{Thoriz}{numeric scalar denoting the time point for which a prediction of the survival status is of interest; \code{Thoriz} must be later than \code{Tstart} and either \code{Dt} or \code{Thoriz} must be specified. If \code{Thoriz} is \code{NULL} is set equal to \code{Tstart + Dt}.}
  
  \item{Dt}{numeric scalar denoting the length of the time interval of prediction; either \code{Dt} or \code{Thoriz} must be specified.}
  
  \item{df_ns}{the degrees of freedom for the natural cubic spline of the cloglog transformation of the predicted
    probabilities used in the Cox model that assesses calibration. The default is 3 unless there are less than 25 events
    in the interval (\code{Tstart}, \code{Thoriz}] in which case it is 2.}
    
  \item{plot}{logical; should a plot be produced. If \code{FALSE}, a list is returned with the observed and predicted probabilities in the calibration plot.}

  \item{add_CI}{logical; should 0.95 pointwise confidence intervals be added around the calibration line.}
  
  \item{col_CI}{character; the color of the shaded area representing the 0.95 pointwise confidence intervals around the calibration line.}
  
  \item{add_density}{logical; should the kernal density estimation of the predicted probabilities be superimposed in the calibration plot.}
  
  \item{col, lwd, lty, col_dens, xlab, ylab, main}{graphical parameters.}
    
  \item{\dots}{additional arguments passed to \code{predict.jm()}.}
}
\value{
If \code{plot = TRUE} (the default) a calibration plot is returned. Otherwise, a list with the following components is returned:
  \item{observed}{numeric vector with observed probabilities.}
  \item{predicted}{grid of probabilities used as x-axis in the calibration plot.}
  \item{pi_u_t}{numeric vector with predicted probabilities.}
  \item{obs_pi_u_t}{transformed predicted probabilities using the Cox model with the cloglog transformation.}
  \item{low}{lower bound of the 95\% confidence interval for the observed probabilities.}
  \item{upp}{upper bound of the 95\% confidence interval for the observed probabilities.}
For the \code{calibration_metrics()} function, a numeric vector with the following components is returned:
  \item{ICI}{numeric scalar denoting the mean absolute difference between the observed and predicted probabilities.}
  \item{E50}{numeric scalar denoting the median absolute difference between the observed and predicted probabilities.}
  \item{E90}{numeric scalar denoting the 90\% percentile of absolute differences between the observed and predicted probabilities.}
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
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
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

jointFit <- jm(CoxFit, fm1, time_var = "year", n_chains = 1L)

# We use landmark time 4 and prediction horizon 3:
calibration_plot(jointFit, newdata = pbc2, Tstart = 4, Dt = 3)

# Now the calibration metrics:
calibration_metrics(jointFit, pbc2, Tstart = 4, Thoriz = 7)
}
}
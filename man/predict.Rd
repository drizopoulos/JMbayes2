\name{Predictions}

\alias{predict.jm}
\alias{plot.predict_jm}

\title{
Predictions from Joint Models
}

\description{
Predict method for object of class \code{"jm"}.
}

\usage{

\method{predict}{jm}(object, newdata = NULL, newdata2 = NULL, times = NULL,
  process = c("longitudinal", "event"),
  type_pred = c("response", "link"),
  type = c("subject_specific", "mean_subject"),
  level = 0.95, return_newdata = FALSE,
  n_samples = 200L, n_mcmc = 55L, cores = NULL, seed = 123L,
  \dots)

\method{plot}{predict_jm}(x, x2 = NULL, subject = 1, outcomes = 1,
  fun_long = NULL, fun_event = NULL, CI = TRUE, xlab = "Follow-up Time",
  ylab_long = NULL, ylab_event = "Cumulative Risk", lwd_long = 2,
  lwd_event = 2, col_line_long = "blue", col_line_event = "red",
  fill_CI_long = "#0000FF80", fill_CI_event = "#FF000080",
  cex_xlab = 1, cex_ylab_long = 1, cex_ylab_event = 1, cex_axis = 1,
  pos_ylab_long = c(0.1, 2, 0.08), \dots)
}

\arguments{
\item{object}{an object inheriting from class \code{"jm"}.}

\item{newdata, newdata2}{data frames}

\item{times}{a numeric vector of future times to calculate predictions.}

\item{process}{for which process to calculation predictions, for the longitudinal outcomes or the event times.}

\item{type_pred}{type of predictions; options are \code{"response"} using the inverse link function in GLMMs, and \code{"link"} that correspond to the linear predictor.}

\item{type}{level of predictions; only relevant when \code{type_pred = "longitudinal"}. Option \code{type = "subject_specific"} combines the fixed- and random-effects parts, whereas \code{type = "mean_subject"} uses only the fixed effects.}

\item{level}{the level of the credible interval.}

\item{return_newdata}{logical; should \code{predict()} return the predictions as extra columns in \code{newdata} and \code{newdata2}.}

\item{n_samples}{the number of samples to use from the original MCMC sample of \code{object}.}

\item{n_mcmc}{the number of Metropolis-Hastings iterations for sampling the random effects per iteration of \code{n_samples}; only the last iteration is retained.}

\item{cores}{how many number of cores to use. If there more than 20 subjects in \code{newdata}, parallel computing is invoked with four cores by default.}

\item{seed}{an integer denoting the seed.}

\item{x, x2}{objects returned by \code{predict.jm()} with argument \code{return_data} set to \code{TRUE}.}

\item{subject}{when multiple subjects are included in the data.frames \code{x} and \code{x2}, it selects which one to plot. Only a single subject can be plotted each time.}

\item{outcomes}{when multiple longitudinal outcomes are included in the data.frames \code{x} and \code{x2}, it selects which ones to plot. A maximum of three outcomes can be plotted each time.}

\item{fun_long, fun_event}{function to apply to the predictions for the longitudinal and event outcomes, respectively. When multiple longitudinal outcomes are plotted, \code{fun_long} can be a list of functions; see examples below.}

\item{CI}{logical; should credible interval areas be plotted.}

\item{xlab, ylab_long, ylab_event}{characture strings or a chracter vector for \code{ylab_long} when multiple longitudinal outcomes are considered with the labels for the horizontal axis, and the two vertical axes.}

\item{lwd_long, lwd_event, col_line_long, col_line_event,
  fill_CI_long, fill_CI_event, cex_xlab, cex_ylab_long, cex_ylab_event,
  cex_axis}{graphical parameters.}

\item{pos_ylab_long}{controls the position of the y-axis labels when multiple longitudinal outcomes are plotted.}

\item{\dots}{extra aguments; currently none is used.}
}

\details{
xxx
}

\value{
Method \code{predict()} returns a list or a data.frame (if \code{return_newdata} was set to \code{TRUE}) with the predictions.

Method \code{plot()} produces figures of the predictions from a single subject.
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{jm}}
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
fm3 <- mixed_model(ascites ~ year * sex, data = pbc2,
                   random = ~ year | id, family = binomial())

jointFit <- jm(CoxFit, list(fm1, fm2, fm3), time_var = "year", n_chains = 2L)

# we select the subject for whom we want to calculate predictions
# we use measurements up to follow-up year 3; we also set that the patients
# were alive up to this time point
t0 <- 3
ND <- pbc2[pbc2$id \%in\% c(2, 25), ]
ND <- ND[ND$year < t0, ]
ND$status2 <- 0
ND$years <- t0

# predictions for the longitudinal outcomes using newdata
predLong1 <- predict(jointFit, newdata = ND, return_newdata = TRUE)

# predictions for the longitudinal outcomes at future time points
# from year 3 to 10
predLong2 <- predict(jointFit, newdata = ND,
                     times = seq(t0, 10, length.out = 51),
                     return_newdata = TRUE)

# predictions for the event outcome at future time points
# from year 3 to 10
predSurv <- predict(jointFit, newdata = ND, process = "event",
                    times = seq(t0, 10, length.out = 51),
                    return_newdata = TRUE)

plot(predLong1)
# for subject 25, outcomes in reverse order
plot(predLong2, outcomes = 3:1, subject = 25)

# prediction for the event outcome
plot(predSurv)

# combined into one plot, the first longitudinal outcome and cumulative risk
plot(predLong2, predSurv, outcomes = 1)

# the first two longitudinal outcomes
plot(predLong1, predSurv, outcomes = 1:2)

# all three longitudinal outcomes, we display survival probabilities instead
# of cumulative risk, and we transform serum bilirubin to the original scale
plot(predLong2, predSurv, outcomes = 1:3, fun_event = function (x) 1 - x,
     fun_long = list(exp, identity, identity),
     ylab_event = "Survival Probabilities",
     ylab_long = c("Serum Bilirubin", "Prothrombin", "Ascites"),
     pos_ylab_long = c(1.9, 1.9, 0.08))
}
}

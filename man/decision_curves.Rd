\name{decision_curves}
\alias{decision_curves}
\alias{decision_curves.default}
\alias{decision_curves.jm}
\alias{decision_curves.coxph}
\alias{plot.decision_curves}

\title{
(Un)weighted Decision Curves for Dynamic Predictions under Left-Truncated and Right-Censored Data
}
\description{
Computes estimates of the net benefit or (weighted) expected cost over a grid of decision thresholds at a specified horizon time (\code{Thoriz}), using longitudinal information available up to a given starting time point (\code{Tstart}). The function supports joint models, Cox proportional hazards models and predictions obtained from any method. It provides inverse probability weighting and model-based estimation methods. For the weighted expected cost, the beta family of weight functions has been implemented. The \code{plot.decision_curves()} method displays the decision curves, supporting multiple models or decision strategies in the same graphic.
}
\usage{
decision_curves(object, \dots)

\method{decision_curves}{default}(object, Time, event, LT_Time = NULL, Tstart, 
                decision_threshold_range = seq(0 + 1e-06, 1 - 1e-06, length.out = 51),
                type = c("EC", "weightedEC", "NB"),
                shape1 = NULL, shape2 = NULL, mod = NULL, interv = NULL, mass = 0.95, 
                Thoriz = NULL, Dt = NULL, type_weights = c("model-based", "IPCW"),
                preds_MB = NULL, opt_interv = c(1, 500), \dots)

\method{decision_curves}{jm}(object, newdata, Tstart, 
                decision_threshold_range = seq(0 + 1e-06, 1 - 1e-06, length.out = 51),
                type = c("EC", "weightedEC", "NB"),
                shape1 = NULL, shape2 = NULL, mod = NULL, interv = NULL, mass = 0.95, 
                Thoriz = NULL, Dt = NULL, 
                type_weights = c("model-based", "IPCW"), opt_interv = c(1, 500), \dots)
      
\method{decision_curves}{coxph}(object, newdata, Tstart, 
                decision_threshold_range = seq(0 + 1e-06, 1 - 1e-06, length.out = 51),
                type = c("EC", "weightedEC", "NB"),
                shape1 = NULL, shape2 = NULL, mod = NULL, interv = NULL, mass = 0.95, 
                Thoriz = NULL, Dt = NULL, type_weights = c("model-based", "IPCW"),
                opt_interv = c(1, 500), \dots)
            
\method{plot}{decision_curves}(x, ...,
     col = c(rep(1, length(objects)), "yellow3", "red"),  
     lty = c(seq_along(objects), 1, 1), 
     lwd = rep(1, length(objects) + 2),
     xlab = "Decision threshold",
     ylab = "Expected cost",
     main = "",
     ylim = NULL, xlim = NULL,
     legend_loc = "topright", legend_cex = 0.8,
     legend_names = c(paste("Model", seq_along(objects)), "Treat all", "Treat none"),
     weight_function = TRUE,
     weight_function_txt = TRUE,
     weight_function_col = adjustcolor("black", alpha.f = 0.6))
}

\arguments{
  \item{object}{an object inheriting from class \code{jm} or \code{coxph}, or a numeric vector containing the dynamic predictions to be evaluated.}

  \item{newdata}{a data.frame that contains the longitudinal and covariate information for the subjects for which prediction of survival probabilities is required. If a \code{jm} object is used, the names of the variables in this data.frame must be the same as in the data.frames that were used to fit the linear mixed effects and the event process model that were supplied as the two first argument of \code{\link{jm}}. If a \code{coxph} is used it only applies to the event process. If the default method is used, the argument is not needed.}
  
  \item{Tstart}{numeric scalar denoting the time point up to which longitudinal information is to be used to derive predictions.}
  
  \item{decision_threshold_range}{a numeric vector containing the grid of decision threshdolds to be used in the decision curves. Default is \code{seq(0 + 1e-06, 1 - 1e-06, length.out = 51)}}
  
  \item{Time, event, LT_Time}{numeric vectors that must be specified when using the default method to compute the metric by directly using a vector of predictions, observed times (\code{Time}), the event indicator (\code{event}) and (possibly) left-truncated times (\code{LT_Time}).}
  
  \item{type}{character string denoting the type of decision curve. Options are expected cost ("EC"; default), weighted expected cost ("weightedEC"; including IWEC estimates) and net benefit ("NB").}
  
  \item{shape1, shape2}{numeric scalars greater or equal to 0 denoting the shape parameteres for the beta family of weight functions (see below for details).}
  
  \item{mod, mass, interv}{two numeric scalars between 0 and 1 and a numeric vector with two scalars between 0 and 1, denoting the mode of the weight function and the mass (default is 95\%) to be allocated within the interval defined by the two values in the numeric vector. These arguments can be used to compute the shape parameters in the weighted decision curves. See below for details.}

  \item{Thoriz}{numeric scalar denoting the time point for which a prediction of the survival status is of interest; \code{Thoriz} must be later than \code{Tstart} and either \code{Dt} or \code{Thoriz} must be specified. If \code{Thoriz} is \code{NULL} is set equal to \code{Tstart + Dt}.}
  
  \item{Dt}{numeric scalar denoting the length of the time interval of prediction; either \code{Dt} or \code{Thoriz} must be specified.}
  
  \item{type_weights}{character string denoting the type of weights to use to account for censorting. Options are model-based (default) and inverse probability weighting (using the Kaplan-Meier estimate of the censoring distribution, and adjusting for left truncation if necessary).}
  
  \item{preds_MB}{a numeric vector containing the model-based weights to be used for the individuals censored within (\code{Tstart}, \code{Thoriz}]. It must be specified when using the default method and the \code{model-based} estimator. The vector must have length equal to the number of individuals considered and must include non-missing elements only for the censored individuals.}

  \item{opt_interv}{a numeric vector containing the end-points of the interval to be searched for the root. This argumentis used when \code{mod}, \code{mass} and \code{interv} are provided, and shape parameters are computed (see details below).}
  
  \item{x}{an object inhereting from class \code{decision_curves} used in the \code{plot} method for decision curves.}
  
  \item{\dots}{additional arguments passed to \code{predict.jm()}. For the \code{plot.decision_curve()} method, it is used to pass more objects of \code{decision_curve} class when multiple decision curves need to be displayed.}
  
  \item{col, lwd, lty, xlab, ylab, xlim, ylim, main}{graphical parameters.}
  
  \item{legend_loc, legend_cex, legend_names}{graphical parameters for the legend displayed in the decision curve.}

  \item{weight_function, weight_function_txt, weight_function_col}{graphical parameters for the weight function displayed in the decision curve when types \code{EC} or \code{weightedEC} are used. \code{weight_function} is a logical; should the weight function be displayed (default is \code{TRUE}). \code{weight_function_txt} is a logical; should the text in the weight function be displayed (default is \code{TRUE}). \code{weight_function_col} is the grapghical parameter for the color.}
}
\details{
The beta family of weight functions is consdiered for the weighted expectec cost. It is defined by the two-parameter form \eqn{w(c;\text{shape1}, \text{shape2})=c^{\text{shape1}-1}(1-c)^{\text{shape2}-1}/B(\text{shape1}, \text{shape2})}, where \eqn{\text{shape1}, \text{shape2}\geq 0} are shape parameters and \eqn{B(\text{shape1},\text{shape2})=1} if \eqn{\text{shape1}=\text{shape2}=0}, and otherwise denotes the beta function.

When \code{mod}, \code{mass} and \code{interv} are provided we compute the shape parameters \code{shape1} and \code{shape2} by solving the system consisting of the two equations: \eqn{\text{mod}=\frac{\text{shape1} - 1}{\text{shape1} + \text{shape2} - 2}} and \eqn{\text{Pr}(X \in \text{interv}) = \text{mass}}, where \eqn{X} follows a beta distribution with shape parameters \code{shape1} and \code{shape2}. To solve this system, the \code{uniroot()} function is internally used. The user can change the \code{interval} argument of the \code{uniroot()} function using the \code{opt_interv} argument.
}
\value{
A list of class \code{decision_curves} with components:
  \item{metric_c}{a list of three components, each containing a numeric vector representing the net benefit or (weighted) expected cost for the prediction model, along with the default strategies of treating all patients or no patients, respectively.}
  \item{newdata}{a data.frame containing the exact dataset used to compute the evaluated dynamic predictions.}
  \item{nr}{a numeric scalar denoting the number of subjects at risk at time \code{Tstart}.}
  \item{nint}{a numeric scalar denoting the number of subjects with an event in (\code{Tstart}, \code{Thoriz}].}
  \item{ncens}{a numeric scalar denoting the number of subjects with a censored time in (\code{Tstart}, \code{Thoriz}].}
  \item{Tstart}{a copy of the \code{Tstart} argument.}
  \item{Thoriz}{a copy of the \code{Thoriz} argument.}
  \item{type}{a copy of the \code{type} argument.}
  \item{shape1, shape2}{a copy of the arguments \code{shape1} and \code{shape2}.}
  \item{decision_threshold_range}{a copy of the \code{decision_threshold_range} argument.}
  \item{iwecs}{a numeric vector containing the IWEC estimates, in case \code{type = "weightedEC"}.}
  \item{type_weights}{a copy of the \code{type_weights} argument.}
  \item{classObject}{the class of \code{object}.}
  \item{nameObject}{the name of \code{object}.}
}
\references{
Vickers, Andrew J and Elkin, Elena B (2006). Decision curve analysis: a novel method for 
evaluating prediction models. \emph{Medical Decision Making} \bold{26}, 565--574.
}
\author{
Arnau Garcia-Fernandez \email{a.garciafernandez@erasmusmc.nl}, Pedro Miranda-Afonso \email{p.mirandaafonso@erasmusmc.nl} and Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}


\seealso{
\code{\link{tvIWEC}}, \code{\link{tvNB}}, \code{\link{tvEC}}, \code{\link{predict}}, \code{\link{jm}}
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

# Classic decision curve using net benefit.
# We use landmark time 4 and prediction horizon 3.
dc_NB <- decision_curves(jointFit, newdata = pbc2, Tstart = 4, Dt = 3, type = "NB")
plot(dc_NB)

# Considering the Cox model, the IPW approach for estimation and 
# shape parameters 12 and 38, we display the decision curve for the expected
# cost, with the weight function at the top part of the plot
dc_EC <- decision_curves(CoxFit, newdata = pbc2.id, Tstart = 2, Thoriz = 6, 
                         shape1 = 12, shape2 = 38, type_weights = "IPCW")
plot(dc_EC, legend_names = c("Cox model", "Treat all strategy", "Treat none strategy"))

# We now display the weighted decision curve. The weight function is such that
# the mode is 0.2, the interval of interest for the decision thresholds is 
# (0.1, 0.3) and the mass to be allocated within the interval is 0.85
dc_wEC <- decision_curves(jointFit, newdata = pbc2, Tstart = 4, Dt = 3, 
                          mod = 0.2, interv = c(0.1, 0.3), mass = 0.8, 
                          type = "weightedEC")
plot(dc_wEC, 
     legend_names = c("Joint model", "Treat all", "Treat none"),
     col = c("black", "lightblue3", "brown2"),
     lty = c(1, 2, 2))
}
}
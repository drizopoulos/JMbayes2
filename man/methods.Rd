\name{jm Methods}
\alias{methods.jm}
\alias{coef}
\alias{coef.jm}
\alias{fixef}
\alias{fixef.jm}
\alias{ranef}
\alias{ranef.jm}
\alias{terms}
\alias{terms.jm}
\alias{model.matrix}
\alias{model.matrix.jm}
\alias{model.frame}
\alias{model.frame.jm}
\alias{family}
\alias{family.jm}
\alias{get_links}
\alias{get_links.jm}
\alias{compare_jm}
\alias{predict}
\alias{predict.jm}

\title{
Various Methods for Standard Generics
}

\description{
Methods for object of class \code{"jm"} for standard generic functions.
}

\usage{

coef(object, \dots)

\method{coef}{jm}(object, \dots)

fixef(object, \dots)

\method{fixef}{jm}(object, outcome = Inf, \dots)

ranef(object, \dots)

\method{ranef}{jm}(object, outcome = Inf, post_vars = FALSE, \dots)

terms(x, \dots)

\method{terms}{jm}(x, process = c("longitudinal", "event"),
                      type = c("fixed", "random"), \dots)

model.frame(formula, \dots)

\method{model.frame}{jm}(formula, process = c("longitudinal", "event"),
                            type = c("fixed", "random"), \dots)

model.matrix(object, \dots)

\method{model.matrix}{jm}(object, \dots)

family(object, \dots)

\method{family}{jm}(object, \dots)

compare_jm(\dots, type = c("marginal", "conditional"),
  order = c("WAIC", "DIC", "LPML", "none"))

predict(object, \dots)

\method{predict}{jm}(object, newdata = NULL, newdata2 = NULL, times = NULL,
  process = c("longitudinal", "event"),
  type_pred = c("response", "link"),
  type = c("subject_specific", "mean_subject"),
  level = 0.95, return_newdata = FALSE,
  n_samples = 200L, n_mcmc = 25L, cores = NULL, seed = 123L,
  \dots)
}

\arguments{
\item{object, x, formula}{object inheriting from class \code{"jm"}.}
\item{outcome}{the index of the linear mixed submodel to extract the estimated fixed effects. If greater
than the total number of submodels, extracts from all of them.}
\item{post_vars}{logical; if \code{TRUE}, returns the variance of the posterior distribution.}
\item{process}{which submodel(s) to extract the terms: \cr
\itemize{
    \item if \code{"longitudinal"}, the linear mixed model(s), or
    \item if \code{"event"}, the survival model.
}
}
\item{type}{in \code{terms()} and \code{model.frame()}, which effects to select in the longitudinal process: \cr
\itemize{
    \item if \code{"fixed"}, the fixed-effects, or
    \item if \code{"random"}, the random-efects.
}
in \code{compare_jm()}, which log-likelihood function use to calculate the criteria: \cr
\itemize{
    \item if \code{"marginal"}, the marginal log-likelihood, or
    \item if \code{"conditional"}, the conditional log-likelihood.
}
in \code{predict()} what type of predictions to calculate.
}
\item{newdata, newdata2}{data frames}
\item{times}{numeric vector of times to predict.}
\item{type_pred}{type of prediction}
\item{level}{level of credible interval.}
\item{return_newdata}{should \code{predict()} return the predictions as extra columns in \code{newdata} and \code{newdata2}.}
\item{n_samples}{number of samples.}
\item{n_mcmc}{number of samples.}
\item{cores}{number of cores for parallel computing.}
\item{seed}{the seed.}
\item{\dots}{further arguments; currently, none is used. \cr in \code{compare_jm()}, a series of \code{jm} objects.}
\item{order}{which criteria use to sort the models in the output.}
}

\details{
\describe{
\item{\code{coef()}}{Extracts estimated fixed effects for the event process from a fitted joint model.}
\item{\code{fixef()}}{Extracts estimated fixed effects for the longitudinal processes from a fitted joint model.}
\item{\code{ranef()}}{Extracts estimated random effects from a fitted joint model.}
\item{\code{terms()}}{Extracts the terms object(s) from a fitted joint model.}
\item{\code{model.frame()}}{Creates the model frame from a fitted joint model.}
\item{\code{model.matrix()}}{Creates the design matrices for linear mixed submodels from a fitted joint model.}
\item{\code{family()}}{Extracts the error distribution and link function used in the linear mixed submodel(s) from a fitted joint model.}
\item{\code{compare_jm()}}{Compares two or more fitted joint models using the criteria WAIC, DIC, and LPML.}
}
}


\value{
\describe{
\item{\code{coef()}}{a list with the elements: \cr
\itemize{
    \item \code{gammas}: estimated baseline fixed effects, and
    \item \code{association}: estimated association parameters.
}
}
\item{\code{fixef()}}{a numeric vector of the estimated fixed effects for the \code{outcome} selected. If the \code{outcome} is greater than the number of linear mixed submodels, it returns a list of numeric vectors for all outcomes.}
\item{\code{ranef()}}{a numeric matrix with rows denoting the individuals and columns the random effects. If \code{postVar = TRUE}, the numeric matrix has the extra attribute "postVar".}
\item{\code{terms()}}{if \code{process = "longitudinal"}, a list of the terms object(s) for the linear mixed model(s).
\cr if \code{process = "event"}, the terms object for the survival model.
}
\item{\code{model.frame()}}{if \code{process = "longitudinal"}, a list of the model frames used in the linear mixed model(s). \cr if \code{process = "event"}, the model frame used in the survival model.
}
\item{\code{model.matrix()}}{a list of the design matrix(ces) for the linear mixed submodel(s).}
\item{\code{family()}}{a list of \code{family} objects.}
\item{\code{compare_jm()}}{a list with the elements: \cr
\itemize{
    \item \code{table}: a table with the criteria calculated for each joint model, and
    \item \code{type}: the log-likelihood function used to calculate the criteria.
}
}
}
}


\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{jm}}
}

\examples{
\donttest{
# linear mixed model fits
fit_lme1 <- lme(log(serBilir) ~ year:sex + age,
                random = ~ year | id, data = pbc2)

fit_lme2 <- lme(prothrombin ~ sex,
                random = ~ year | id, data = pbc2)

# cox model fit
fit_cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

# joint model fit
fit_jm <- jm(fit_cox, list(fit_lme1, fit_lme2), time_var = "year",
    n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)

# coef(): fixed effects for the event process
coef(fit_jm)

# fixef(): fixed effects for the first linear mixed submodel
fixef(fit_jm, outcome = 1)

# ranef(): random effects from all linear mixed submodels
head(ranef(fit_jm))

# terms(): random effects terms for the first linear mixed submodel
terms(fit_jm, process = "longitudinal", type = "random")[[1]]

# mode.frame(): model frame for the fixed effects in the second
# linear mixed submodel
head(model.frame(fit_jm, process = "longitudinal", type = "fixed")[[2]])

# model.matrix(): fixed effects design matrix for the first linear
# mixed submodel
head(model.matrix(fit_jm)[[1]])

# family(): family objects from both linear mixed submodels
family(fit_jm)

# compare_jm(): compare two fitted joint models
fit_lme1b <- lme(log(serBilir) ~ 1,
                  random = ~ year | id, data = pbc2)

fit_jm2 <- jm(fit_cox, list(fit_lme1b, fit_lme2), time_var = "year",
    n_chains = 1L, n_iter = 11000L, n_burnin = 1000L)

compare_jm(fit_jm, fit_jm2)
}
}

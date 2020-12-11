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
}

\arguments{
\item{object, x, formula}{object inheriting from class \code{"jm"}.}
\item{outcome}{the linear mixed submodel to extract the estimated fixed effects. If greater
than the total number of submodels, extracts from all of them.}
\item{post_vars}{logical; if \code{TRUE} returns the variance of the posterior distribution.}
\item{process}{which submodel(s) to extract the terms:
\tabular{ll}{
\code{"longitudinal"} \tab linear mixed model(s), \cr
\code{"survival"} \tab survival model.}
}
\item{type}{in \code{compare_jm()}, which effects to select in the longitudinal process:
\tabular{ll}{
\code{"fixed"} \tab fixed-effects, \cr
\code{"random"} \tab random-efects.} \cr
in \code{terms()}, which log-likelihood function use to calculate the criteria:
\tabular{ll}{
\code{"marginal"} \tab marginal log-likelihood, \cr
\code{"conditional"} \tab conditional log-likelihood.} \cr
}
\item{\dots}{further arguments; currently none is used. \cr\cr in \code{compare_jm()}, a series of \code{jm} objects.}
\item{order}{which criteria use to sort the models in the output.}
}

\value{
\tabular{ll}{

\code{coef()} \tab a list with the elements: 
  \tabular{ll}{
  \code{gammas} \tab estimated baseline fixed effects. \cr
  \code{association} \tab estimated association parameters.} \cr

\code{fixef()} \tab a numeric vector of the estimated fixed effects for the \code{outcome} selected. If \code{outcome} is greater than the number of linear mixed submodels, returns a list of numeric vectors for all outcomes. \cr

\code{ranef()} \tab a numeric matrix with rows denoting the individuals and columns the random effects. If \code{postVar = TRUE}, the numeric matrix has an extra attribute "postVar". \cr

\code{terms()} \tab 
\tabular{ll}{ 
if \code{process = "longitudinal"}, a list of the terms object(s) for the linear mixed model(s). \cr if \code{process = "event"}, the terms object for the survival model.} \cr

\code{model.frame()} \tab 
\tabular{ll}{
if \code{process = "longitudinal"}, a list of the model frames used in the linear mixed model(s). \cr If \code{process = "event"}, the model frame used in the survival model.} \cr

\code{model.matrix()} \tab a list of the design matrices for the linear mixed submodels. \cr

\code{family()} \tab a list of \code{family} objects. \cr

\code{compare_jm()} \tab an object of class \code{compare_jm}. This is a list with the elements: 
\tabular{ll}{
\code{table} \tab a table with the criteria calculted to each joint model. \cr
\code{type} \tab the log-likelihood function used to calculate the criteria.   
} \cr
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
fit_jm <- jm(fit_cox, list(fit_lme1, fit_lme2), time_var = "year")

# coef(): fixed effects for the event process 
coef(fit_jm)

# fixef(): fixed effects for the first linear mixed submodels
fixef(fit_jm, outcome = 1)

# ranef(): random effects from all linear mixed submodels
ranef(fit_jm)

# terms(): fixed effects terms in the linear mixed models
terms(fit_jm, process = "longitudinal", type = "random")

# mode.frame(): model frame for the fixed effects terms in the linear mixed 
# submodels
model.frame(fit_jm, process = "longitudinal", type = "fixed")

# model.matrix(): linear mixed models' design matrices
model.matrix(fit_jm)

# family(): family objects from all linear mixed submodels
family(fit_jm)

# compare_jm(): compare two fitted joint models
fit_lme1_b <- lme(log(serBilir) ~ year:sex,
                 random = ~ year | id, data = pbc2)

fit_jm_b <- jm(fit_cox, list(fit_lme1_b, fit_lme2), time_var = "year")

compare_jm(fit_jm, fit_jm_b)
}
}

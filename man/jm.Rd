\name{jm}
\alias{jm}
\alias{value}
\alias{slope}
\alias{area}
\alias{tv}

\title{
Joint Models for Longitudinal and Time-to-Event Data
}

\description{
Fits multivariate joint models for longitudinal and time-to-event data.
}

\usage{
jm(Surv_object, Mixed_objects, time_var, functional_forms = NULL,
  data_Surv = NULL, id_var = NULL, priors = NULL, control = NULL, \dots)

value(x)
slope(x)
area(x)

tv(x, knots = NULL, ord = 2L)
}

\arguments{
\item{Surv_object}{an object: \cr 
\itemize{
    \item of class 'coxph' fitted by function \code{coxph()} from package \bold{survival} or
    \item of class 'survreg' fitted by function \code{survreg()} from package \bold{survival}
    }}
\item{Mixed_objects}{an object: \cr 
\itemize{
    \item of class 'lme' fitted by function \code{lme()} from package \bold{nlme} or
    \item of class 'MixMod' fitted by function \code{mixed_model()} from package \bold{GLMMadaptive}
    }}
\item{time_var}{a character string indicating the time variable in the mixed effects model.}
\item{functional_forms}{a list of formulas. Each formula corresponds to one longitudinal outcome and specifies the association structure between that outcome and the survival submodel as well as interaction terms between the components of the longitudinal outcome and the survival submodel. See \bold{Examples}.}
\item{data_Surv}{the data.frame used to fit the Cox/AFT survival submodel.}
\item{id_var}{a character string indicating the id variable in the survival submodel.}
\item{priors}{.}
\item{control}{.}
\item{x}{.}
\item{knots, ord}{.}
\item{\dots}{arguments passed to \code{control}.}
}

\details{
XXX
}

\value{
XXX
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{methods.jm}}
}


\examples{
x <- 5L
}

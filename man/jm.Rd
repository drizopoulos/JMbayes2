\name{jm}
\alias{jm}
\alias{value}
\alias{slope}
\alias{area}

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
}

\arguments{
\item{Surv_object}{.}
\item{Mixed_objects}{.}
\item{time_var}{.}
\item{functional_forms}{.}
\item{data_Surv}{.}
\item{id_var}{.}
\item{priors}{.}
\item{control}{.}
\item{x}{.}
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

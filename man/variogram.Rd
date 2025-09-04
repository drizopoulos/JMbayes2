\name{variogram}
\alias{variogram}

\title{ Variogram for Longitudinal Data}

\description{
  It computes the semi-variogram for longitudinal data.
}

\usage{
variogram(y, times, id, compute_var = TRUE)
}

\arguments{
  \item{y}{a numeric vector of longitudinal responses.}
  \item{times}{a numeric vector of times at which the longitudinal responses were collected.}
  \item{id}{a numeric vector of a factor of subject id numbers.}
  \item{compute_var}{logical; if \code{TRUE} the variance is also computed.}
}

\value{
  A list with two components, i.e., \code{svar}, a two-column matrix with the time lags and the variogram
  values, and \code{sigma2} the total variance.
}

\author{
    Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\examples{
ind <- aids$patient %in% c(1, 2)
yy <- aids$CD4[ind]
tt <- aids$obstime[ind]
ids <- aids$patient[ind]
variogram(yy, tt, ids)
}

\keyword{methods}

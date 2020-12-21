\name{crLong}
\alias{crLong}

\title{ Transform Competing Risks Data in Long Format }

\description{
  In a competing risks setting this function expands the data frame with a
  single row per subject to the a data frame in long format in which each
  subject has as many rows as the number of competing events.
}

\usage{
crLong(data, statusVar, censLevel,
    nameStrata = "strata", nameStatus = "status2")
}

\arguments{
  \item{data}{the data frame containing the competing risk data with a single
      row per subject.}
  \item{statusVar}{a character string denoting the name of the variable in
          \code{data} that identifies the status variable which equals 1 if the
          subject had any of the competing events and 0 otherwise.}
  \item{censLevel}{a character string or a scalar denoting the censoring level
      in the \code{statusVar} variable of \code{data}.}
  \item{nameStrata}{a character string denoting the variable that will be added
    in the long version of \code{data} denoting the various causes of event.}
  \item{nameStatus}{a character string denoting the variable that will be added
    in the long version of \code{data} denoting if the subject experience any
    of the competing events.}
}

\value{
  A data frame in the long format with multiple rows per subject.
}

\references{
Rizopoulos, D. (2012) \emph{Joint Models for Longitudinal and Time-to-Event Data: with 
Applications in R}. Boca Raton: Chapman and Hall/CRC.

Putter, H., Fiocco, M., and Geskus, R. (2007). Tutorial in biostatistics:
Competing risks and multi-state models. \emph{Statistics in Medicine} \bold{26},
2389--2430.
}

\author{
    Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\examples{
head(crLong(pbc2.id, "status", "alive"))
}

\keyword{methods}

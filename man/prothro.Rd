\name{prothro}
\alias{prothro}
\alias{prothros}

\docType{data}

\title{Prednisone versus Placebo in Liver Cirrhosis Patients}

\description{
  A randomized trial on 488 liver cirrhosis patients
}

\format{
  Two data frames with the following variable.
  \describe{
    \item{\code{id}}{patients identifier; in total there are 467 patients.}
    \item{\code{pro}}{prothrobin measurements.}
    \item{\code{time}}{for data frame \code{prothro} the time points at which the prothrobin measurements were taken; 
        for data frame \code{prothros} the time to death or censoring.}
    \item{\code{death}}{a numeric vector with 0 denoting censoring and 1 death.}
    \item{\code{treat}}{randomized treatment; a factor with levels "placebo" and "prednisone".}
  }
}

\source{
  \url{http://www.gllamm.org/books/readme.html#14.6}, 
}

\references{

Andersen, P. K., Borgan, O., Gill, R. D. and Keiding, N. (1993). 
\emph{Statistical Models Based on Counting Processes}. New York: Springer. 

}

\keyword{datasets}

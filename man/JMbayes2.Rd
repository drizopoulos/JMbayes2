\name{JMbayes2}
\alias{JMbayes2-package}
\alias{JMbayes2}
\docType{package}

\title{
Extended Joint Models for Longitudinal and Time-to-Event Data
}

\description{
Fit joint models for longitudinal and time-to-event data under the Bayesian approach. Multiple longitudinal outcomes of mixed type (continuous/categorical) and multiple event times (competing risks and multi-state processes) are accommodated.
}

\details{
\tabular{ll}{
Package: \tab JMbayes2\cr
Type: \tab Package\cr
Version: \tab 0.1-7\cr
Date: \tab 2021-07-05\cr
License: \tab GPL (>=3)\cr
}

This package fits joint models for longitudinal and time-to-event data. It can accommodate multiple longitudinal outcomes of different type (e.g., continuous, dichotomous, ordinal, counts), and assuming different distributions, i.e., Gaussian, Student's-t, Gamma, Beta, unit Lindley, censored Normal, Binomial, Poisson, Negative Binomial, and Beta-Binomial. For the event time process, right, left and interval censored data can be handled, while competing risks and multi-sate processes are also covered.

\pkg{JMbayes2} fits joint models using Markov chain Monte Carlo algorithms implemented in C++. The package also offers several utility functions that can extract useful information from
fitted joint models. The most important of those are included in the
\bold{See also} Section below.
}

\author{
Dimitris Rizopoulos, Grigorios Papageorgiou, Pedro Miranda Afonso

Maintainer: Dimitris Rizopoulos <d.rizopoulos@erasmusmc.nl>
}

\references{
Rizopoulos, D. (2012). Joint Models for Longitudinal and Time-to-Event Data With Applications in R. Boca Raton: Chapman & Hall/CRC.
}

\keyword{package}
\keyword{multivariate}

\seealso{
    \code{\link{jm}},
    \code{\link{methods.jm}},
    \code{\link{coda_methods.jm}}
}


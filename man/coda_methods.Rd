\name{jm coda Methods}
\alias{coda_methods.jm}
\alias{traceplot}
\alias{traceplot.jm}
\alias{gelman_diag}
\alias{gelman_diag.jm}
\alias{densplot}
\alias{densplot.jm}
\alias{cumuplot}
\alias{cumuplot.jm}

\title{
Various Methods for Functions from the coda package
}

\description{
Methods for object of class \code{"jm"} for standard generic functions.
}

\usage{

traceplot(object, \dots)

\method{traceplot}{jm}(object,
  parm = c("all", "betas", "sigmas", "D", "bs_gammas",
           "tau_bs_gammas", "gammas", "alphas"), \dots)

gelman_diag(object, \dots)

\method{gelman_diag}{jm}(object,
  parm = c("all", "betas", "sigmas", "D", "bs_gammas",
           "tau_bs_gammas", "gammas", "alphas"), \dots)

densplot(object, \dots)

\method{densplot}{jm}(object,
  parm = c("all", "betas", "sigmas", "D", "bs_gammas",
           "tau_bs_gammas", "gammas", "alphas"), \dots)

cumuplot(object, \dots)

\method{cumuplot}{jm}(object,
  parm = c("all", "betas", "sigmas", "D", "bs_gammas",
           "tau_bs_gammas", "gammas", "alphas"), \dots)

}

\arguments{
\item{object}{objects inheriting from class \code{"jm"}.}
\item{parm}{.}
\item{\dots}{further arguments passed to the corresponding function from the
  the \bold{coda} package.}
}

\details{
XXX
}


\value{
XXX.
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{jm}}
}

\examples{
\donttest{
x <- rnorm(5L)
}
}

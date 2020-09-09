\name{jm coda Methods}
\alias{coda_methods.jm}
\alias{traceplot}
\alias{traceplot.jm}
\alias{ggtraceplot}
\alias{ggtraceplot.jm}
\alias{gelman_diag}
\alias{gelman_diag.jm}
\alias{densplot}
\alias{densplot.jm}
\alias{ggdensityplot}
\alias{ggdensityplot.jm}
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

ggtraceplot(object, \dots)

\method{ggtraceplot}{jm}(object,
  parm = c("all", "betas", "sigmas", "D", "bs_gammas",
           "tau_bs_gammas", "gammas", "alphas"),
  size = 1, alpha = 0.8,
  chaincols = c('standard', 'catalog', 'metro',
                'pastel', 'beach', 'moonlight', 'goo', 'sunset'),
  gridrow = 3, gridcol = 1, grid = FALSE, \dots)

gelman_diag(object, \dots)

\method{gelman_diag}{jm}(object,
  parm = c("all", "betas", "sigmas", "D", "bs_gammas",
           "tau_bs_gammas", "gammas", "alphas"), \dots)

densplot(object, \dots)

\method{densplot}{jm}(object,
  parm = c("all", "betas", "sigmas", "D", "bs_gammas",
           "tau_bs_gammas", "gammas", "alphas"), \dots)

ggdensityplot(object, \dots)

\method{ggdensityplot}{jm}(object,
  parm = c("all", "betas", "sigmas", "D", "bs_gammas",
           "tau_bs_gammas", "gammas", "alphas"),
  size = 1, alpha = 0.6,
  chaincols = c('standard', 'catalog', 'metro', 'pastel',
                'beach', 'moonlight', 'goo', 'sunset'),
  gridrow = 3, gridcol = 1, grid = FALSE, \dots)

cumuplot(object, \dots)

\method{cumuplot}{jm}(object,
  parm = c("all", "betas", "sigmas", "D", "bs_gammas",
           "tau_bs_gammas", "gammas", "alphas"), \dots)

}

\arguments{
\item{object}{objects inheriting from class \code{"jm"}.}
\item{parm}{.}
\item{size}{.}
\item{alpha}{.}
\item{chaincols}{.}
\item{gridrow}{.}
\item{gridcol}{.}
\item{grid}{.}
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

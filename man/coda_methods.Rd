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
\item{object}{an object inheriting from class \code{"jm"}.}
\item{parm}{a character string specifying which parameters of the joint model to plot. Possible values are:  Defaults to \code{'all'}, \code{'betas'}, \code{'alphas'}, \code{'sigmas'}, \code{'D'}, \code{'bs_gammas'}, \code{'tau_bs_gammas'}, \code{'gammas'}.\cr}
\item{size}{the width of the traceplot line in mm. Defaults to 1.}
\item{alpha}{the opacity level of the traceplot line. Defaults to 0.8.}
\item{theme}{a character string specifying the color theme to be used. Possible options are \code{'standard'}, \code{'catalog'}, \code{'metro'}, \code{'pastel'}, \code{'beach'}, \code{'moonlight'}, \code{'goo'}, \code{'sunset'}. \cr}
\item{grid}{logical; defaults to \code{FALSE}. If \code{TRUE} the plots are returned in grids split over multiple pages. For more details see the documentation for \code{\link[gridExtra:arrangeGrob]{gridExtra::marrangeGrob()}}. \cr}
\item{gridrows}{number of rows per page for the grid. Only relevant when using \code{grid = TRUE}. Defaults to 3.}
\item{gridcols}{number of columns per page for the grid. Only relevant when using \code{grid = TRUE}. Defaults to 1.}
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

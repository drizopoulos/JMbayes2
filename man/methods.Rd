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

\method{fixef}{jm}(object, outcome = 1, \dots)

ranef(object, \dots)

\method{ranef}{jm}(object, outcome = Inf, post_vars = FALSE, \dots)

terms(x, \dots)

\method{terms}{jm}(x, process = c("longitudinal", "event"),
                      type = c("fixed", "random"), \dots)

model.frame(formula, \dots)

\method{model.frame}{jm}(formula, process = c("longitudinal", "event"),
                            type = c("fixed", "random"), \dots)

model.matrix(object, \dots)

\method{model.matrix}{jm}(object, process = c("longitudinal", "event"),
                             type = c("fixed", "random"), \dots)

family(object, \dots)

\method{family}{jm}(object, \dots)
}

\arguments{
\item{object, x, formula}{objects inheriting from class \code{"jm"}.}
\item{outcome}{.}
\item{post_vars}{.}
\item{process}{.}
\item{type}{.}
\item{\dots}{further arguments; currently none is used.}
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

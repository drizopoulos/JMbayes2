\name{lme}
\alias{lme}
\alias{coxph}
\alias{mixed_model}

\title{Slice-aware model fitting generics}

\description{
These functions are S3 generics exported by \pkg{JMbayes2} to enable method dispatch
on the class of the \code{data} argument. When \code{data} is a \code{"sliced_data"}
object (as produced by \code{\link{slicer}}), the corresponding slice-wise fitting
methods are used.
}

\usage{
lme(fixed, data, ...)
coxph(formula, data, ...)
mixed_model(fixed, data, ...)
}

\arguments{
\item{fixed}{For \code{lme()} and \code{mixed_model()}, a model formula for the fixed effects.}

\item{formula}{For \code{coxph()}, a survival model formula.}

\item{data}{A data frame (default methods) or a \code{"sliced_data"} object (slice-wise methods).}

\item{...}{Further arguments passed to the underlying model-fitting functions.
For slice-wise methods, additional arguments include \code{parallel_out} and \code{cores}
to control parallel execution across slices.}
}

\details{
When \code{data} is a regular data frame, the default methods call
\code{nlme::lme()}, \code{survival::coxph()}, and \code{GLMMadaptive::mixed_model()},
respectively.

When \code{data} is a \code{"sliced_data"} object, the corresponding \code{*.sliced_data}
methods fit the requested model independently within each slice and return a list of
fits (one per slice).
}

\note{
JMbayes2 exports S3 generics \code{lme()}, \code{coxph()}, and \code{mixed_model()}
to enable dispatch on \code{"sliced_data"}. When JMbayes2 is attached, these names
mask \code{nlme::lme}, \code{survival::coxph}, and \code{GLMMadaptive::mixed_model}.
Use the \code{pkg::fun} form to call the original functions.
}

\seealso{
\code{\link{slicer}}, \code{nlme::lme}, \code{survival::coxph}, \code{GLMMadaptive::mixed_model}.
}

\examples{
\dontrun{

slc <- slicer(n_slices = 2, id_var = "id", data_long = pbc2, data_surv = pbc2.id)
n_cores <- max(parallel::detectCores() - 1L, 1L)

lme_fit <- lme(fixed  = log(serBilir) ~ year * sex,
               data   = slc$long,
               random = ~ year | id,
               cores  = n_cores)

cox_fit <- coxph(formula = Surv(years, status2) ~ sex,
                 data    = slc$surv,
                 cores   = n_cores)

mxm_fit <- mixed_model(fixed  = ascites ~ year + sex,
                       data   = slc$long,
                       random = ~ year | id,
                       family = binomial(),
                       cores  = n_cores)
}
}

\keyword{methods}

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
Various Methods for Functions from the \emph{coda} Package
}

\description{
Methods for an object of class \code{"jm"} for diagnostic functions.
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
  theme = c('standard', 'catalog', 'metro',
                'pastel', 'beach', 'moonlight', 'goo', 'sunset', 'custom'),
  grid = FALSE, gridrows = 3, gridcols = 1, custom_theme = NULL, \dots)

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
  theme = c('standard', 'catalog', 'metro', 'pastel',
                'beach', 'moonlight', 'goo', 'sunset', 'custom'),
  grid = FALSE, gridrows = 3, gridcols = 1, custom_theme = NULL, \dots)

cumuplot(object, \dots)

\method{cumuplot}{jm}(object,
  parm = c("all", "betas", "sigmas", "D", "bs_gammas",
           "tau_bs_gammas", "gammas", "alphas"), \dots)

}

\arguments{
\item{object}{an object inheriting from class \code{"jm"}.}
\item{parm}{a character string specifying which parameters of the joint model to plot. Possible options are \code{'all'}, \code{'betas'}, \code{'alphas'}, \code{'sigmas'}, \code{'D'}, \code{'bs_gammas'}, \code{'tau_bs_gammas'}, or \code{'gammas'}.\cr}
\item{size}{the width of the traceplot line in mm. Defaults to 1.}
\item{alpha}{the opacity level of the traceplot line. Defaults to 0.8.}
\item{theme}{a character string specifying the color theme to be used. Possible options are \code{'standard'}, \code{'catalog'}, \code{'metro'}, \code{'pastel'}, \code{'beach'}, \code{'moonlight'}, \code{'goo'}, or \code{'sunset'}. Note that this option supports fitted objects with up to three chains. If you have an object with more chains then the colors are either automatically chosen, or you can specify your own custom color theme via argument \code{custom_theme}. \cr}
\item{grid}{logical; defaults to \code{FALSE}. If \code{TRUE}, the plots are returned in grids split over multiple pages. For more details see the documentation for \code{\link[gridExtra:arrangeGrob]{gridExtra::marrangeGrob()}}. \cr}
\item{gridrows}{number of rows per page for the grid. Only relevant when using \code{grid = TRUE}. Defaults to 3.}
\item{gridcols}{number of columns per page for the grid. Only relevant when using \code{grid = TRUE}. Defaults to 1.}
\item{custom_theme}{A named character vector with elements equal to the number of chains (\code{n_chains}). The name of each element should be the number corresponding to the respective chain. Defaults to \code{NULL}.}
\item{\dots}{further arguments passed to the corresponding function from the \bold{coda} package.}
}

\value{
\describe{
\item{\code{traceplot()}}{Plots the evolution of the estimated parameter vs. iterations in a fitted joint model.}

\item{\code{ggtraceplot()}}{Plots the evolution of the estimated parameter vs. iterations in a fitted joint model using \pkg{ggplot2}.}

\item{\code{gelman_diag()}}{Calculates the potential scale reduction factor for the estimated parameters in a fitted joint model, together with the upper confidence limits.}

\item{\code{densplot()}}{Plots the density estimate for the estimated parameters in a fitted joint model.}

\item{\code{ggdensityplot()}}{Plots the evolution of the estimated parameter vs. iterations in a fitted joint model using \pkg{ggplot2}.}

\item{\code{cumuplot()}}{Plots the evolution of the sample quantiles vs. iterations in a fitted joint model.}
}
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{jm}}
}

\examples{
\donttest{
# linear mixed model fits
fit_lme1 <- lme(log(serBilir) ~ year:sex + age,
                random = ~ year | id, data = pbc2)

fit_lme2 <- lme(prothrombin ~ sex,
                random = ~ year | id, data = pbc2)

# cox model fit
fit_cox <- coxph(Surv(years, status2) ~ age, data = pbc2.id)

# joint model fit
fit_jm <- jm(fit_cox, list(fit_lme1, fit_lme2), time_var = "year", n_chains = 1L)

# trace plot for the fixed effects in the linear mixed submodels
traceplot(fit_jm, parm = "betas")

# density plot for the fixed effects in the linear mixed submodels
densplot(fit_jm, parm = "betas")

# cumulative quantile plot for the fixed effects in the linear mixed submodels
cumuplot(fit_jm, parm = "betas")

# trace plot for the fixed effects in the linear mixed submodels
ggtraceplot(fit_jm, parm = "betas")
ggtraceplot(fit_jm, parm = "betas", grid = TRUE)
ggtraceplot(fit_jm, parm = "betas", custom_theme = c('1' = 'black'))

# trace plot for the fixed effects in the linear mixed submodels
ggdensityplot(fit_jm, parm = "betas")
ggdensityplot(fit_jm, parm = "betas", grid = TRUE)
ggdensityplot(fit_jm, parm = "betas", custom_theme = c('1' = 'black'))
}
}


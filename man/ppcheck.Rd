\name{ppcheck}
\alias{ppcheck}

\title{ Posterio Predictive Checks for Joint Models}

\description{
  It computes various posterior predictive checks for joint models.
}

\usage{
ppcheck(object, nsim = 40L, newdata = NULL, seed = 123L,
    process = c("longitudinal", "event"),
    type = c("ecdf", "variogram", "variance_function"),
    CI_ecdf = c("binomial", "Dvoretzky–Kiefer–Wolfowitz"),
    outcomes = Inf, percentiles = c(0.025, 0.975),
    random_effects = c("posterior_means", "mcmc", "prior"),
    Fforms_fun = NULL, ylim = NULL, \dots)
}

\arguments{
  \item{object}{an object inheriting from class \code{"jm"}.}
  \item{nsim}{a numeric scalar denoting the number of replicated data.}
  \item{newdata}{a data.frame based on which to simulate replicated data.
        Relevant for cross-validated posterior predictive checks.}
  \item{seed}{the seed to use.}
  \item{process}{a character string indicating for which process to do the checks.}
  \item{type}{a character string indicating the type of checks. The default \code{ecdf} compares the
        empirical distribution function of the replicated data with the one of the observed data. For
        continous longitudinal outcomes, the semi-variogram can also be used.}
  \item{CI_ecdf}{character string indicating the type of confidence interval for the
        emprical distribution function.}
  \item{outcomes}{a numeric vector of indices indicating for which longitudinal outcomes to do the checks.
        The deault value \code{Inf} implies checking all longitudinal outcomes.}
  \item{percentiles}{a numeric vector of length two indicating the percentiles of the observed data to use when
        depicting the checks.}
  \item{random_effects}{a character string indicating what to do with the random effects.}
  \item{Fforms_fun}{a function that calculates the functional forms.}
  \item{ylim}{the \code{ylim} argument for the plot.}
  \item{\dots}{extra argument passed to the \code{simualte} method.}
}

\value{
  a figure with the checks.
}

\author{
    Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\examples{
\donttest{
# Cox model for the composite event death or transplantation
pbc2.id$status2 <- as.numeric(pbc2.id$status != 'alive')
pbc2$status2 <- as.numeric(pbc2$status != 'alive')
CoxFit <- coxph(Surv(years, status2) ~ sex, data = pbc2.id)

# a linear mixed model for log serum bilirubin
fm1 <- lme(log(serBilir) ~ ns(year, 3) * sex, data = pbc2,
           random = list(id = pdDiag(~ ns(year, 3))))

# the joint model
jointFit <- jm(CoxFit, fm1, time_var = "year", save_random_effects = TRUE)

ppcheck(jointFit)

FF <- function (t, betas, bi, data) {
    sex <- as.numeric(data$sex == "female")
    NS <- ns(t, k = c(0.9911, 3.9863), B = c(0, 14.10579))
    X <- cbind(1, NS, sex, NS * sex)
    Z <- cbind(1, NS)
    eta <- c(X \%*\% betas[[1]]) + rowSums(Z * bi)
    cbind(eta)
}

ppcheck(jointFit, process = "event", Fforms_fun = FF)

}
}

\keyword{methods}

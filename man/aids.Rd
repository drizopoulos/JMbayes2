\name{aids}
\alias{aids}
\alias{aids.id}

\docType{data}

\title{Didanosine versus Zalcitabine in HIV Patients}

\description{
  A randomized clinical trial in which both longitudinal and survival data were collected to compare the efficacy and 
  safety of two antiretroviral drugs in treating patients who had failed or were intolerant of zidovudine (AZT) therapy.
}

\format{
  A data frame with 1408 observations on the following 9 variables.
  \describe{
    \item{\code{patient}}{patients identifier; in total there are 467 patients.}
    \item{\code{Time}}{the time to death or censoring.}
    \item{\code{death}}{a numeric vector with 0 denoting censoring and 1 death.}
    \item{\code{CD4}}{the CD4 cells count.}
    \item{\code{obstime}}{the time points at which the CD4 cells count was recorded.}
    \item{\code{drug}}{a factor with levels \code{ddC} denoting zalcitabine and \code{ddI} denoting didanosine.}
    \item{\code{gender}}{a factor with levels \code{female} and \code{male}.}
    \item{\code{prevOI}}{a factor with levels \code{AIDS} denoting previous opportunistic infection (AIDS 
        diagnosis) at study entry, and \code{noAIDS} denoting no previous infection.}
    \item{\code{AZT}}{a factor with levels \code{intolerance} and \code{failure} denoting AZT intolerance and 
        AZT failure, respectively.}
  }
}

\references{
Goldman, A., Carlin, B., Crane, L., Launer, C., Korvick, J., Deyton, L. and Abrams, D. (1996) Response of CD4+ and
clinical consequences to treatment using ddI or ddC in patients with advanced HIV infection. \emph{Journal of Acquired
Immune Deficiency Syndromes and Human Retrovirology} \bold{11}, 161--169.

Guo, X. and Carlin, B. (2004) Separate and joint modeling of longitudinal and event time data using standard
computer packages. \emph{The American Statistician} \bold{58}, 16--24.
}

\note{
The data frame \code{aids.id} contains the first CD4 cell count measurement for each patient. This data frame is used to 
fit the survival model.
}

\keyword{datasets}

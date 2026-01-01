\name{consensus}
\alias{consensus}

\title{ Consensus Combination of Posterior Draws from Joint Model Fits }

\description{
  This function combines posterior draws from multiple joint model fits into a 
  single set of draws.
}

\usage{
consensus(object, parm, method = c("union", "equal_weight", "var_weight"),
          seed = 123L)
}

\arguments{
  \item{object}{an object of class \code{"sliced_jm"}.}
  \item{parm}{a character vector with the names of parameters to combine. These 
              must correspond to elements in \code{fit$mcmc[[parm]]} for each fit
              (e.g., \code{"gammas"}, \code{"alphas"}, \code{"betas1"}, \code{"betas2"}).}
  \item{method}{the consensus method used to combine draws.
  \describe{
    \item{\code{"union"}}{concatenate draws across slices (no averaging).}
    \item{\code{"equal_weight"}}{compute an iteration-wise simple average across slices.}
    \item{\code{"var_weight"}}{compute an iteration-wise weighted average across slices using 
                               inverse-variance weights.}
  }}
  \item{seed}{an integer seed used for the within-slice random permutation step 
              (used by \code{"equal_weight"} and \code{"var_weight"}).}
}

\value{
An object of class \code{"consensus_jm"} with components:
\describe{
  \item{\code{method}}{the selected method.}
  \item{\code{parm}}{the requested parameter(s) block(s).}
  \item{\code{n_splits}}{number of fits combined.}
  \item{\code{draws}}{named list of combined draw matrices (one per parameter block).}
  \item{\code{weights}}{for \code{"equal_weight"} and \code{"var_weight"}, a named list of weight matrices.}
  \item{\code{n_draws}}{number of combined draws per parameter block.}
  \item{\code{summary}}{named list of summary matrices with \code{Mean}, \code{StDev}, \code{2.5\%}, \code{97.5\%}, and \code{P}).}
  \item{\code{seed}}{seed used.}
}
}

\author{
    Pedro Miranda-Afonso \email{p.mirandaafonso@erasmusmc.nl}
}

\keyword{methods}

\name{slicer}
\alias{slicer}

\title{ Split Longitudinal and Survival Data into Subject-level Samples }

\description{
  This function partitions a longitudinal dataset and a survival dataset at the
  subject level, returning two lists of sliced datasets for model fitting in 
  parallel or series.
}

\usage{
slicer(n_slices, id_var, data_long, data_surv, seed = 123L)
}

\arguments{
  \item{n_slices}{an integer scalar giving the number of data slices (subsamples) to create.}
  \item{id_var}{a character scalar with the name of the subject identifier 
                variable in both \code{data_long} and \code{data_surv}.}
  \item{data_long}{a data frame containing the longitudinal measurements.}
  \item{data_surv}{a data frame containing the survival (time-to-event) information.}
  \item{seed}{an integer seed used to randomize the assignment of subject IDs to slices.}
}

\value{
A list with two components:
\describe{
  \item{\code{long}}{a list of length \code{n_slices} with class \code{"sliced_data"}. 
                     Each element is a data frame containing the longitudinal 
                     rows for the subjects assigned to that slice.}
  \item{\code{surv}}{a list of length \code{n_slices} with class \code{"sliced_data"}. 
                     Each element a data frame containing the survival rows for 
                     the subjects assigned to that slice.}
}
}

\examples{
data(pbc2, package = "JMbayes2")
data(pbc2.id, package = "JMbayes2")

pbc2_slc <- slicer(n_slices = 2, id_var = "id", data_long = pbc2, 
                   data_surv = pbc2.id, seed = 123L)
length(pbc2_slc$long) # 2
length(pbc2_slc$surv) # 2
}

\author{
    Pedro Miranda-Afonso \email{p.mirandaafonso@erasmusmc.nl}
}

\keyword{methods}

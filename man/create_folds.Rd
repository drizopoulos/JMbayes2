\name{create_folds}
\alias{create_folds}

\title{
Creating Training and Testing Datasets 
}
\description{
This is an auxiliar function that splits the data into folds for internal validation or Bootstrap methods.
}
\usage{
create_folds(data, V = 5, id_var = "id",
    method = c("CV", "Bootstrap"), strata = NULL, seed = 123L)
}
\arguments{
  \item{data}{the data.frame to split in folds.}
  
  \item{V}{numeric scalar denoting the number of folds for cross-validation or the number of sample for the Bootstrap methods.}
  
  \item{id_var}{character string denoting the name of the subject id variable in \code{data}.}
  
  \item{method}{character string indicating which method to use to create the training and testing datasets in
  \code{create_folds()}. The default is V-fold cross-validation. For the \code{Bootstrap} option, \code{V} samples
  with replacement from the original dataset are proruced as training data. The testing data contains the subjects
  that were not selected in the respective Bootstrap sample.}
  
  \item{strata}{character vector with the names of stratifying variables.}
  
  \item{seed}{integer denoting the seed.}
}
\value{
A list with the following components:
  \item{training}{list containing training datasets.}
  \item{testing}{list containing testing datasets.}
}
\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}
\seealso{
\code{\link{predict}}, \code{\link{jm}}
}
\examples{
\donttest{
CVdats <- create_folds(pbc2, V = 5, id_var = "id")
CVdats
}
}


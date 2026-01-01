\name{rc_setup}
\alias{rc_setup}

\title{ Combine Recurring and Terminal Event Data in Long Format }

\description{
  This function combines two data frames, the recurring-event and terminal-event/competing-risks 
  datasets, into one. Each subject has as many rows in the new data frame as the 
  number of recurrent risk periods plus one for each terminal event/competing risk.
}

\usage{
rc_setup(rc_data, trm_data,
    idVar = "id", statusVar = "status",
    startVar = "start", stopVar = "stop",
    trm_censLevel,
    nameStrata = "strata", nameStatus = "status")
}

\arguments{
  \item{rc_data}{the data frame containing the recurring-event data with 
      multiple rows per subject.}
  \item{trm_data}{the data frame containing the terminal-event/competing-risks 
      data with a single row per subject.}
  \item{idVar}{a character scalar denoting the name of the variable in
      \code{rc_data} and \code{trm_data} that identifies the subject/group.}
  \item{statusVar}{a character string denoting the name of the variable in
      \code{rc_data} and \code{trm_data} that identifies the status variable. 
      In \code{rc_data} equals 1 if the subject had an event and 0 otherwise.
      In \code{trm_data} equals to the event or censoring level.}
   \item{startVar}{a character string denoting the name of the variable in
      \code{rc_data} that identifies the starting time for the risk interval.}
    \item{stopVar}{a character string denoting the name of the variable in
      \code{rc_data} and \code{trm_data} that identifies the event or censoring 
      time.}
    \item{trm_censLevel}{a character string or a scalar denoting the censoring 
    level in the statusVar variable of \code{trm_data}.}
    \item{nameStrata}{a character string denoting the variable that will be added
      in the long version of \code{data} denoting the various causes of event.}
    \item{nameStatus}{a character string denoting the variable that will be added
      in the long version of \code{data} denoting if the subject had an event.}
}

\value{
  A data frame in the long format with multiple rows per subject.
}

\author{
    Pedro Miranda-Afonso \email{p.mirandaafonso@erasmusmc.nl}
}

\keyword{methods}

\name{rc_setup}
\alias{rc_setup}

\title{ Combine Recurring and Terminal Event Data in Long Format }

\description{
  This function combines two data frames, the recurring event and terminal event 
  datasets, into one. Each subject has as many rows in the new data frame as the 
  number of recurrent risk periods plus one for the terminal event.
}

\usage{
rc_setup(rc_data, trm_data,
    rc_idVar = "id", rc_statusVar = "status",
    rc_startVar = "start", rc_stopVar = "stop",
    trm_idVar = "id", trm_statusVar = "status",
    trm_stopVar = "stop",
    nameStrata = "strata", nameStatus = "status")
}

\arguments{
  \item{rc_data}{the data frame containing the recurring event data with 
      multiple rows per subject.}
  \item{trm_data}{the data frame containing the terminal event data with 
      a single row per subject.}
  \item{rc_idVar}{a character string denoting the name of the variable in
      \code{rc_data} that identifies the subject.}
  \item{rc_statusVar}{a character string denoting the name of the variable in
      \code{rc_data} that identifies the status variable which equals 1 if the
      subject had an event and 0 otherwise.}
   \item{rc_startVar}{a character string denoting the name of the variable in
      \code{rc_data} that identifies the starting time for the interval.}
    \item{rc_stopVar}{a character string denoting the name of the variable in
      \code{rc_data} that identifies the stopping time for the interval.}
    \item{trm_idVar}{a character string denoting the name of the variable in
      \code{trm_data} that identifies the subject.}
    \item{trm_statusVar}{a character string denoting the name of the variable in
      \code{trm_data} that identifies the status variable which equals 1 if the
      subject had the event and 0 otherwise.}
    \item{trm_stopVar}{a character string denoting the name of the variable in
      \code{trm_data} that identifies the follow up time.}
    \item{nameStrata}{a character string denoting the variable that will be added
      in the long version of \code{data} denoting the various causes of event.}
    \item{nameStatus}{a character string denoting the variable that will be added
      in the long version of \code{data} denoting if the subject had an event.}
}

\value{
  A data frame in the long format with multiple rows per subject.
}

\author{
    Pedro Miranda Afonso \email{p.mirandaafonso@erasmusmc.nl}
}

\keyword{methods}

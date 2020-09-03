\name{pbc2}
\alias{pbc2}
\alias{pbc2.id}


\docType{data}

\title{Mayo Clinic Primary Biliary Cirrhosis Data}

\description{
Followup of 312 randomised patients with primary biliary cirrhosis, a rare autoimmune liver disease, at Mayo Clinic.
}

\format{
  A data frame with 1945 observations on the following 20 variables.
  \describe{
    \item{\code{id}}{patients identifier; in total there are 312 patients.}
    \item{\code{years}}{number of years between registration and the earlier of death, transplantion, or study 
        analysis time.}
    \item{\code{status}}{a factor with levels \code{alive}, \code{transplanted} and \code{dead}.}
    \item{\code{drug}}{a factor with levels \code{placebo} and \code{D-penicil}.}
    \item{\code{age}}{at registration in years.}
    \item{\code{sex}}{a factor with levels \code{male} and \code{female}.}
    \item{\code{year}}{number of years between enrollment and this visit date, remaining values on the line of 
        data refer to this visit.}
    \item{\code{ascites}}{a factor with levels \code{No} and \code{Yes}.}
    \item{\code{hepatomegaly}}{a factor with levels \code{No} and \code{Yes}.}
    \item{\code{spiders}}{a factor with levels \code{No} and \code{Yes}.}
    \item{\code{edema}}{a factor with levels \code{No edema} (i.e., no edema and no diuretic therapy for edema),
        \code{edema no diuretics} (i.e., edema present without diuretics, or edema resolved by diuretics), and 
        \code{edema despite diuretics} (i.e., edema despite diuretic therapy).}
    \item{\code{serBilir}}{serum bilirubin in mg/dl.}
    \item{\code{serChol}}{serum cholesterol in mg/dl.}
    \item{\code{albumin}}{albumin in gm/dl.}
    \item{\code{alkaline}}{alkaline phosphatase in U/liter.}
    \item{\code{SGOT}}{SGOT in U/ml.}
    \item{\code{platelets}}{platelets per cubic ml / 1000.}
    \item{\code{prothrombin}}{prothrombin time in seconds.}
    \item{\code{histologic}}{histologic stage of disease.}
    \item{\code{status2}}{a numeric vector with the value 1 denoting if the patient was dead,
        and 0 if the patient was alive or transplanted.}
  }
}

\references{
Fleming, T. and Harrington, D. (1991) \emph{Counting Processes and Survival Analysis}. Wiley, New York.

Therneau, T. and Grambsch, P. (2000) \emph{Modeling Survival Data: Extending the Cox Model}. Springer-Verlag, New York.
}

\note{
The data frame \code{pbc2.id} contains the first measurement for each patient. This data frame is used to 
fit the survival model. 
}

\keyword{datasets}

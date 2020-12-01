\name{jm}
\alias{jm}
\alias{value}
\alias{slope}
\alias{area}
\alias{tv}

\title{
Joint Models for Longitudinal and Time-to-Event Data
}

\description{
Fits multivariate joint models for longitudinal and time-to-event data.
}

\usage{
jm(Surv_object, Mixed_objects, time_var, functional_forms = NULL,
  data_Surv = NULL, id_var = NULL, priors = NULL, control = NULL, \dots)

value(x)
slope(x)
area(x)

tv(x, knots = NULL, ord = 2L)
}

\arguments{
\item{Surv_object}{an object: \cr 
\itemize{
    \item of class 'coxph' fitted by function \code{coxph()} from package \bold{survival} or
    \item of class 'survreg' fitted by function \code{survreg()} from package \bold{survival}
    }}
\item{Mixed_objects}{an object: \cr 
\itemize{
    \item of class 'lme' fitted by function \code{lme()} from package \bold{nlme} or
    \item of class 'MixMod' fitted by function \code{mixed_model()} from package \bold{GLMMadaptive}
    }}
\item{time_var}{a character string indicating the time variable in the mixed effects model.}
\item{functional_forms}{a list of formulas. Each formula corresponds to one longitudinal outcome and specifies the association structure between that outcome and the survival submodel as well as interaction terms between the components of the longitudinal outcome and the survival submodel. See \bold{Examples}.}
\item{data_Surv}{the data.frame used to fit the Cox/AFT survival submodel.}
\item{id_var}{a character string indicating the id variable in the survival submodel.}
\item{priors}{a named list of user-specified prior parameters:
  \describe{
        \item{mean_betas_HC}{}
        \item{Tau_betas_HC}{}
        \item{mean_betas_nHC}{}
        \item{Tau_betas_nHC}{}
        \item{mean_bs_gammas}{the prior mean vector of the normal prior for the B-splines 
            coefficients used to approximate the baseline hazard.}
        \item{Tau_bs_gammas}{the prior precision matrix of the normal prior for the B-splines 
            coefficients used to approximate the baseline hazard.}
        \item{A_tau_bs_gammas}{the prior shape parameter of the Gamma prior for the
            precision parameter of the penalty term for the B-splines coefficients for
            the baseline hazard.}
        \item{B_tau_bs_gammas}{the prior rate parameter of the Gamma prior for the
            precision parameter of the penalty term for the B-splines coefficients for
            the baseline hazard.}
        \item{rank_Tau_bs_gammas}{}
        \item{mean_gammas}{the prior mean vector of the normal prior for the regression 
            coefficients of baseline covariates.}
        \item{Tau_gammas}{the prior precision matrix of the normal prior for the regression 
            coefficients of baseline covariates.}
        \item{penalty_gammas}{a character string with value: 'none', 'ridge', 'horseshoe' indicating whether the coefficients of the baseline covariates included in the survival submodel should not be shrinked, shrinked using ridge prior or shrinked using horseshoe prior respectively.}
        \item{A_lambda_gammas}{}
        \item{B_lambda_gammas}{}
        \item{A_tau_gammas}{}
        \item{B_tau_gammas}{}
        \item{A_nu_gammas}{}
        \item{B_nu_gammas}{}
        \item{A_xi_gammas}{}
        \item{B_xi_gammas}{}
        \item{mean_alphas}{the prior mean vector of the normal prior for the association
            parameters.}
        \item{Tau_alphas}{the prior mean vector of the normal prior for the association
            parameters.}
        \item{penalty_alphas}{a character string with value: 'none', 'ridge', 'horseshoe' indicating whether the coefficients association parameters should not be shrinked, shrinked using ridge prior or shrinked using horseshoe prior respectively.}
        \item{A_lambda_alphas}{}
        \item{B_lambda_alphas}{}
        \item{A_tau_alphas}{}
        \item{B_tau_alphas}{}
        \item{A_nu_alphas}{}
        \item{B_nu_alphas}{}
        \item{A_xi_alphas}{}
        \item{B_xi_alphas}{}
        \item{prior_D_sds_df}{}
        \item{prior_D_sds_sigma}{}
        \item{prior_D_L_etaLKJ}{}
        \item{prior_sigmas_df}{}
        \item{prior_sigmas_sigma}{}
    }
}
\item{control}{.}
\item{x}{.}
\item{knots, ord}{.}
\item{\dots}{arguments passed to \code{control}.}
}

\details{
XXX
}

\value{
XXX
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{methods.jm}}
}


\examples{
x <- 5L
}

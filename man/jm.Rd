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
    \item of class 'coxph' fitted by function \code{coxph()} from package \bold{survival}, or
    \item of class 'survreg' fitted by function \code{survreg()} from package \bold{survival}.
    }}
\item{Mixed_objects}{a \code{list} of objects or a single object. Objects may be: \cr 
\itemize{
    \item of class 'lme' fitted by function \code{lme()} from package \bold{nlme}, or
    \item of class 'MixMod' fitted by function \code{mixed_model()} from package \bold{GLMMadaptive}.
    }}
\item{time_var}{a \code{character} string indicating the time variable in the mixed-effects model(s).}
\item{functional_forms}{a \code{list} of formulas. Each formula corresponds to one longitudinal outcome and specifies the association structure between that outcome and the survival submodel as well as any interaction terms between the components of the longitudinal outcome and the survival submodel. See \bold{Examples}.}
\item{data_Surv}{the \code{data.frame} used to fit the Cox/AFT survival submodel.}
\item{id_var}{a \code{character} string indicating the id variable in the survival submodel.}
\item{priors}{a named \code{list} of user-specified prior parameters:
  \describe{
        \item{\code{mean_betas_HC}}{the prior mean vector of the normal prior for the regression coefficients of the covariates of the longitudinal model(s), which were hierarchically centered.}
        \item{\code{Tau_betas_HC}}{the prior precision matrix of the normal prior for the regression coefficients of the longitudinal model(s), which were hierarchically centered.}
        \item{\code{mean_betas_nHC}}{a \code{list} of the prior mean vector(s) of the normal prior(s) for the regression coefficients of the covariates of the longitudinal model(s), which were not hierarchically centered.}
        \item{\code{Tau_betas_nHC}}{a \code{list} of the prior precision matrix(ces) of the normal prior(s) for the regression coefficients of the longitudinal model(s), which were not Hierarchically Centered.}
        \item{\code{mean_bs_gammas}}{the prior mean vector of the normal prior for the B-splines 
            coefficients used to approximate the baseline hazard.}
        \item{\code{Tau_bs_gammas}}{the prior precision matrix of the normal prior for the B-splines 
            coefficients used to approximate the baseline hazard.}
        \item{\code{A_tau_bs_gammas}}{the prior shape parameter of the gamma prior for the
            precision parameter of the penalty term for the B-splines coefficients for
            the baseline hazard.}
        \item{\code{B_tau_bs_gammas}}{the prior rate parameter of the gamma prior for the
            precision parameter of the penalty term for the B-splines coefficients for
            the baseline hazard.}
        \item{\code{rank_Tau_bs_gammas}}{the prior rank parameter for the precision matrix of the normal prior for the B-splines coefficients used to approximate the baseline hazard.}
        \item{\code{mean_gammas}}{the prior mean vector of the normal prior for the regression 
            coefficients of baseline covariates.}
        \item{\code{Tau_gammas}}{the prior precision matrix of the normal prior for the regression 
            coefficients of baseline covariates.}
        \item{\code{penalty_gammas}}{a character string with value 'none', 'ridge', or 'horseshoe' indicating whether the coefficients of the baseline covariates included in the survival submodel should not be shrunk, shrank using ridge prior, or shrank using horseshoe prior, respectively.}
        \item{\code{A_lambda_gammas}}{the prior shape parameter of the gamma prior for the
            precision parameter of the local penalty term for the baseline regression coefficients. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{B_lambda_gammas}}{the prior rate parameter of the gamma prior for the
            precision parameter of the local penalty term for the baseline regression coefficients. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{A_tau_gammas}}{the prior shape parameter of the gamma prior for the
            precision parameter of the global penalty term for the baseline regression coefficients. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{B_tau_gammas}}{the prior rate parameter of the gamma prior for the
            precision parameter of the global penalty term for the baseline regression coefficients. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{A_nu_gammas}}{the prior shape parameter of the gamma prior for the variance hyperparameter for the precision parameter of the local penalty term for the baseline regression coefficients. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{B_nu_gammas}}{the prior rate parameter of the gamma prior for the variance hyperparameter for the precision parameter of the local penalty term for the baseline regression coefficients. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{A_xi_gammas}}{the prior shape parameter of the gamma prior for the variance hyperparameter for the precision parameter of the global penalty term for the baseline regression coefficients. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{B_xi_gammas}}{the prior rate parameter of the gamma prior for the variance hyperparameter for the precision parameter of the global penalty term for the baseline regression coefficients. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{mean_alphas}}{the prior mean vector of the normal prior for the association
            parameter(s).}
        \item{\code{Tau_alphas}}{the prior mean vector of the normal prior for the association
            parameter(s).}
        \item{\code{penalty_alphas}}{a character string with value 'none', 'ridge', 'horseshoe' indicating whether the coefficients association parameters should not be shrunk, shrank using ridge prior, or shrank using horseshoe prior, respectively.}
        \item{\code{A_lambda_alphas}}{the prior shape parameter of the gamma prior for the
            precision parameter of the local penalty term for the association parameters. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{B_lambda_alphas}}{the prior rate parameter of the gamma prior for the
            precision parameter of the local penalty term for the association parameters. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{A_tau_alphas}}{the prior shape parameter of the gamma prior for the
            precision parameter of the global penalty term for the association parameters. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{B_tau_alphas}}{the prior rate parameter of the gamma prior for the
            precision parameter of the global penalty term for the association parameters. Only relevant when \code{penalty_gammas = 'ridge'} or \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{A_nu_alphas}}{the prior shape parameter of the gamma prior for the variance hyperparameter for the precision parameter of the local penalty term for the association parameters. Only relevant when \code{penalty_gammas = 'ridge'}, or \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{B_nu_alphas}}{the prior rate parameter of the gamma prior for the variance hyperparameter for the precision parameter of the local penalty term for the association parameters. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{A_xi_alphas}}{the prior shape parameter of the gamma prior for the variance hyperparameter for the precision parameter of the global penalty term for the association parameters. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{B_xi_alphas}}{the prior rate parameter of the gamma prior for the variance hyperparameter for the precision parameter of the global penalty term for the association parameters. Only relevant when \code{penalty_gammas = 'ridge'} or when \code{penalty_gammas = 'horseshoe'}.}
        \item{\code{gamma_prior_D_sds}}{logical; if \code{TRUE}, a gamma prior will be used for the standard deviations of the D matrix (variance-covariance matrix of the random effects); default is \code{TRUE}}
        \item{\code{D_sds_df}}{the prior degrees of freedom parameter for the half-t prior for the standard deviations of the D matrix (variance-covariance matrix of the random effects).}
        \item{\code{D_sds_sigma}}{the prior sigma parameter vector for the half-t prior for the standard deviations of the D matrix (variance-covariance matrix of the random effects).}
        \item{\code{D_sds_shape}}{the prior shape parameter for the gamma prior for the standard deviations of the D matrix (variance-covariance matrix of the random effects).}
        \item{\code{D_sds_mean}}{the prior mean parameter vector for the gamma prior for the standard deviations of the D matrix (variance-covariance matrix of the random effects).}
        \item{\code{D_L_etaLKJ}}{the prior eta parameter for the LKJ prior for the correlation matrix of the random effects.}
        \item{\code{sigmas_df}}{the prior degrees of freedom parameter for the half-t prior for the error term(s).}
        \item{\code{sigmas_sigma}}{the prior sigma parameter for the half-t prior for the error term(s).}
    }
}
\item{control}{a list of control values with components:
    \describe{        
        \item{\code{GK_k}}{the number of quadrature points for the Gauss Kronrod rule; options 15 and 7.}
        \item{\code{Bsplines_degree}}{the degree of the splines in each basis; default quadratic splines.}
        \item{\code{base_hazard_segments}}{the number of segments to split the follow-up period; default is 10.}
        \item{\code{diff}}{the order of the difference used in the penalty matrix for the B-splines for h_0; default is 2.}
        \item{\code{n_chains}}{an integer specifying the number of chains for the MCMC; default is 3.}
        \item{\code{n_burnin}}{an integer specifying the number of burn-in iterations; default is 500.}
        \item{\code{n_iter}}{an integer specifying the number of total iterations per chain; default is 3500.}
        \item{\code{seed}}{the seed used in the sampling procedures; default is 123.}
        \item{\code{MALA}}{a \code{logical}; if TRUE, the MALA algorithm is used when updating the elements of the Cholesky factor of the D matrix; default is \code{FALSE}.}
        \item{\code{save_random_effects}}{a \code{logical}; if TRUE, the full MCMC results of the random effects will be saved and returned with the \code{jm} object; default is \code{FALSE}.}
        \item{\code{cores}}{an integer specifying the number of cores to use; default is the available cores minus one.}
    }
  }
\item{x}{.}
\item{knots, ord}{.}
\item{\dots}{arguments passed to \code{control}.}
}

\details{
XXX
}

\value{
A list of class \code{jm} with components:

\item{mcmc}{a \code{list} of the MCMC samples for each parameter.}
\item{acc_rates}{a \code{list} of the acceptance rates for each parameter.}
\item{logLik}{a \code{matrix} of dimensions [\code{((n_iter - n_burnin)/n_thin)*n_thin}, number of individuals], with element [i, j] being the conditional log-Likelihood value of the \eqn{i^{th}} iteration for the \eqn{j^{th}} individual.}
\item{mlogLik}{a \code{matrix} of dimensions [\code{((n_iter - n_burnin)/n_thin)*n_thin}, number of individuals], with element [i, j] being the marginal log-Likelihood value of the \eqn{i^{th}} iteration for the \eqn{j^{th}} individual.}
\item{running_time}{an object of class \code{proc_time} with the time used to run \code{jm}.}
\item{statistics}{a \code{list} with posterior estimates of the parameters (means, medians, standard deviations, standard errors, effective sample sizes, tail probabilities, upper and lower bounds of credible intervals, etc.).}
\item{fit_stats}{a \code{list} of lists with fit statistics (DIC, pD, LPML, CPO, WAIC) for both conditional and marginal formulations.}
\item{model_data}{a \code{list} of data used to fit the model.}
\item{model_info}{a \code{list} of components of the fit useful to other functions.}
\item{initial_values}{a \code{list} with the initial values of the parameters.}
\item{control}{a copy of the \code{control} values used to fit the model.}
\item{priors}{a copy of the \code{priors} used to fit the model.}
\item{call}{the matched call.}
}

\author{
Dimitris Rizopoulos \email{d.rizopoulos@erasmusmc.nl}
}

\seealso{
\code{\link{methods.jm}}
}


\examples{
\dontrun{
################################################################################

##############################################
# Univariate joint model for serum bilirubin #
# 1 continuous outcome                       #                       
##############################################

# [1] Fit the mixed model using lme(). 
fm1 <- lme(fixed = log(serBilir) ~ year * sex + I(year^2) + 
           age + prothrombin, random =  ~ year | id, data = pbc2)

# [2] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [3] The basic joint model is fitted using a call to jm() i.e.,
joint_model_fit_1 <- jm(fCox1, fm1, time_var = "year")
summary(joint_model_fit_1)
traceplot(joint_model_fit_1)

################################################################################

##########################################################################
# Multivariate joint model for serum bilirubin, hepatomegaly and ascites # 
# 1 continuous outcome, 2 categorical outcomes                           #
##########################################################################

# [1] Fit the mixed-effects models using lme() for continuous 
# outcomes and mixed_model() for categorical outcomes. 
fm1 <- lme(fixed = log(serBilir) ~ year * sex, 
           random = ~ year | id, data = pbc2)

fm2 <- mixed_model(hepatomegaly ~ sex + age + year, data = pbc2, 
                   random = ~ year | id, family = binomial())

fm3 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ year | id, family = binomial())

# [2] Save all the fitted mixed-effects models in a list.
Mixed <- list(fm1, fm2, fm3)

# [3] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model. 
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [4] The joint model is fitted using a call to jm() i.e.,
joint_model_fit_2 <- jm(fCox1, Mixed, time_var = "year")
summary(joint_model_fit_2)
traceplot(joint_model_fit_2)

################################################################################

######################
# Slope & Area Terms # 
######################

# We extend model 'joint_model_fit_2' by including the value and slope term for
# bilirubin, the area term for hepatomegaly (in the log-odds scale), and the
# value and area term for spiders (in the log-odds scale). 
# To include these terms into the model, we specify the 'functional_forms' 
# argument. This should be a list of right side formulas. Each component of the
# list should have as name the name of the corresponding outcome variable. In
# the right side formula we specify the functional form of the association using 
# functions 'value()', 'slope()' and 'area()'.
# Notes: (1) For terms not specified in the 'functional_forms' list, the default
# value functional form is used. 

# [1] Fit the mixed-effects models using lme() for continuous outcomes 
# and mixed_model() for categorical outcomes. 
fm1 <- lme(fixed = log(serBilir) ~ year * sex, random = ~ year | id, data = pbc2)

fm2 <- mixed_model(hepatomegaly ~ sex + age + year, data = pbc2, 
                   random = ~ year | id, family = binomial())

fm3 <- mixed_model(ascites ~ year + age, data = pbc2,
                   random = ~ year | id, family = binomial())

# [2] Save all the fitted mixed-effects models in a list.
Mixed <- list(fm1, fm2, fm3)

# [3] Fit a Cox model, specifying the baseline covariates to be included in the
# joint model.
fCox1 <- coxph(Surv(years, status2) ~ drug + age, data = pbc2.id)

# [4] Specify the list of formulas to be passed to the functional_forms argument 
# of jm().
fForms <- list("log(serBilir)" = ~ value(log(serBilir)) + slope(log(serBilir)),
               "hepatomegaly" = ~ area(hepatomegaly),
               "ascites" = ~ value(ascites) + area(ascites))
               
# [5] The joint model is fitted using a call to jm() and passing the list 
# to the functional_forms argument.
joint_model_fit_2 <- jm(fCox1, Mixed, time_var = "year", 
                        functional_forms = fForms)
summary(joint_model_fit_2)

}
}

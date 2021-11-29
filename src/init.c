#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _JMbayes2_cum_haz(SEXP, SEXP);
extern SEXP _JMbayes2_hSfun(SEXP, SEXP);
extern SEXP _JMbayes2_logLik_jm(SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes2_mcmc_cpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes2_mlogLik_jm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _JMbayes2_simulate_REs(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_JMbayes2_cum_haz",      (DL_FUNC) &_JMbayes2_cum_haz,      2},
    {"_JMbayes2_hSfun",        (DL_FUNC) &_JMbayes2_hSfun,        2},
    {"_JMbayes2_logLik_jm",    (DL_FUNC) &_JMbayes2_logLik_jm,    4},
    {"_JMbayes2_mcmc_cpp",     (DL_FUNC) &_JMbayes2_mcmc_cpp,     5},
    {"_JMbayes2_mlogLik_jm",   (DL_FUNC) &_JMbayes2_mlogLik_jm,   6},
    {"_JMbayes2_simulate_REs", (DL_FUNC) &_JMbayes2_simulate_REs, 3},
    {NULL, NULL, 0}
};

void R_init_JMbayes2(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

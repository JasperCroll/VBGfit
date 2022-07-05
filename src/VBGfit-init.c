/*
# calls for the c code of the VBGfit package
# 15 march 2022
# latest version on https://bitbucket.org/JCroll/vbgfit/
*/
 
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


/* .Call calls */
extern SEXP minLL_growthcurve(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP predict_growthcurve(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"minLL_growthcurve",   (DL_FUNC) &minLL_growthcurve,   5},
    {"predict_growthcurve", (DL_FUNC) &predict_growthcurve, 4},
};

void R_init_popstruct(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

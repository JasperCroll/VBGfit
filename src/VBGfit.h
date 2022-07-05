/*
specific header file for the C implementations of the VBGfit package
15 march 2022
latest version on https://bitbucket.org/JCroll/vbgfit/
*/

#include <R.h>
#include <Rdefines.h>

#ifndef __VBGfit_H__
#define __VBGfit_H__

double predict_length_avg(double lprev, double fRlm, double rB, int schrinking);
double predict_length_var(double lprev, double fRlm, double rB);

double calculate_prob(double len, double len_avg, double len_var);

SEXP predict_growthcurve(SEXP poppars, SEXP globpars, SEXP modeltype, SEXP feedingbounds);

SEXP minLL_growthcurve(SEXP poppars, SEXP globpars, SEXP samples, SEXP modeltype, SEXP feedingbounds);

#endif

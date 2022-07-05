/*
This file contains the C implementation of routines in the VBGfit package
15 march 2022
latest version on https://bitbucket.org/JCroll/vbgfit/
*/
 
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "VBGfit.h"

// calculate average size of cohort in next timestep
double predict_length_avg(double lprev, double fRlm, double rB, int schrinking){
    double lnext;
    if(schrinking == 0 && lprev > fRlm){
        lnext = lprev;
    }
    else{
        lnext = lprev*exp(-rB)+(1-exp(-rB))*fRlm;
    }
    return lnext;
}

// calculate variance in size of cohort in next timestep
double predict_length_var(double lprev, double fRlm, double rB){
    double lnext = lprev*exp(-2*rB)+0.5*rB*(1-exp(-2*rB))*fRlm;
    return lnext;
}

double calculate_prob(double len, double len_avg, double len_var){
    double prob = 1/( sqrt( 2*M_PI * len_var ) ) * exp( -0.5*pow( ( len-len_avg),2 ) / len_var );
          //Rprintf("break = %.20f %f %.20f, %f %f, %f %f \n", prob, 1/( sqrt( 2*M_PI * len_var ) ) , exp( -0.5*pow( ( len-len_avg),2 ) / len_var ), -0.5*pow( ( len-len_avg),2 ), len_var, len, len_avg     );
    if(prob<1E-100) prob = 1E-100;

    return prob;
}


// predict growth curve over the years
SEXP predict_growthcurve(SEXP poppars, SEXP globpars, SEXP modeltype, SEXP feedingbounds){

    int avg = 0;
    int var = 1;

    int model = asInteger(modeltype);

    //split parameters
    int npars = INTEGER(globpars)[0];
    int nyears = INTEGER(globpars)[1];
    int nages = INTEGER(globpars)[3];
    int schrinking = INTEGER(globpars)[6];
    int logscale = INTEGER(globpars)[7];
    int grouppar = INTEGER(globpars)[8];
    int ngroups = INTEGER(globpars)[9];
    

    double *** fRlm = (double***) malloc(2* sizeof(double**));
    double ** len0 = (double**) malloc(2* sizeof(double*));
    double ** year0 = (double**) malloc(2* sizeof(double*));
    double rB;
    
    fRlm[avg] = (double**) malloc(ngroups*sizeof(double*));
    fRlm[var] = (double**) malloc(ngroups*sizeof(double*));
    
    len0[avg] = (double*) malloc(nyears*sizeof(double));
    len0[var] = (double*) malloc(nyears*sizeof(double));
    year0[avg] = (double*) malloc(nages*sizeof(double));
    year0[var] = (double*) malloc(nages*sizeof(double));
    
    for(int i=0; i<ngroups; i++){
        fRlm[avg][i] = (double*) malloc((nyears-1)*sizeof(double));
        fRlm[var][i] = (double*) malloc((nyears-1)*sizeof(double));
    }
    

    // split popparameter values
    if(model == 0 && logscale == 1){
        if(npars != 2*ngroups+2*nyears+2*(nages-1)+1) error("incorect number of parameters for model with constant feeding level (C code)");
        
        for(int i=0;i<nyears;i++){
            if(i<nyears-1){
                for(int j=0; j<ngroups; j++){
                    fRlm[avg][j][i] = exp((double)REAL(poppars)[j]);
                    fRlm[var][j][i] = exp((double)REAL(poppars)[ngroups+j]);
                }
            }
            
            len0[avg][i] = exp((double)REAL(poppars)[2*ngroups+i]);
            len0[var][i] = exp((double)REAL(poppars)[2*ngroups+nyears+i]);
        }
        
        year0[avg][0] = exp((double)REAL(poppars)[2*ngroups]);
        year0[var][0] = exp((double)REAL(poppars)[2*ngroups+nyears]);
        
        for(int i=1;i<nages;i++){
            year0[avg][i] = exp((double)REAL(poppars)[2*ngroups + 2*nyears + i-1]);
            year0[var][i] = exp((double)REAL(poppars)[2*ngroups + 2*nyears + nages-1 + i-1]);
        }
        rB = exp((double)REAL(poppars)[2*ngroups + 2*nyears + 2*(nages-1)]);
    }
    else if(model == 0 && logscale == 0){
        if(npars != 2*ngroups+2*nyears+2*(nages-1)+1) error("incorect number of parameters for model with constant feeding level (C code)");
        
        for(int i=0;i<nyears;i++){
            if(i<nyears-1){
                for(int j=0; j<ngroups; j++){
                    fRlm[avg][j][i] = (double)REAL(poppars)[j];
                    fRlm[var][j][i] = (double)REAL(poppars)[ngroups+j];
                }
            }
            
            len0[avg][i] = (double)REAL(poppars)[2*ngroups+i];
            len0[var][i] = (double)REAL(poppars)[2*ngroups+nyears+i];
        }
        
        year0[avg][0] = (double)REAL(poppars)[2*ngroups];
        year0[var][0] = (double)REAL(poppars)[2*ngroups+nyears];
        
        for(int i=1;i<nages;i++){
            year0[avg][i] = (double)REAL(poppars)[2*ngroups + 2*nyears + i-1];
            year0[var][i] = (double)REAL(poppars)[2*ngroups + 2*nyears + nages-1 + i-1];
        }
        rB = (double)REAL(poppars)[2*ngroups + 2*nyears + 2*(nages-1)];
    }
    else if(model == 1 && logscale == 1){

        if(npars != 2*ngroups*(nyears-1) + 2*nyears + 2*(nages-1)  + 1) error("incorect number of parameters for model with varying feeding level (C code)");

        for(int i=0;i<nyears;i++){
            if(i<nyears-1){
                for(int j=0; j<ngroups; j++){
                    fRlm[avg][j][i] = exp((double)REAL(poppars)[j*(nyears-1)+i]);
                    fRlm[var][j][i] = exp((double)REAL(poppars)[ngroups*(nyears-1)+j*(nyears-1)+i]);
                }
            }
            len0[avg][i] = exp((double)REAL(poppars)[2*ngroups*(nyears-1)+i]);
            len0[var][i] = exp((double)REAL(poppars)[2*ngroups*(nyears-1) + nyears +i]);
        }
        year0[avg][0] = exp((double)REAL(poppars)[2*ngroups*(nyears-1)]);
        year0[var][0] = exp((double)REAL(poppars)[2*ngroups*(nyears-1) + nyears]);
        for(int i=1;i<nages;i++){
            year0[avg][i] = exp((double)REAL(poppars)[2*ngroups*(nyears-1) + 2*nyears + i-1]);
            year0[var][i] = exp((double)REAL(poppars)[2*ngroups*(nyears-1) + 2*nyears + nages-1 + i-1]);
        }
        rB = exp((double)REAL(poppars)[2*ngroups*(nyears-1) + 2*nyears + 2*(nages-1) ]);

    }
    else if(model == 1 && logscale == 0){
        
        if(npars !=  2*ngroups*(nyears-1) + 2*nyears + 2*(nages-1)  + 1) error("incorect number of parameters for model with varying feeding level (C code)");

        for(int i=0;i<nyears;i++){
            if(i<nyears-1){
                for(int j=0; j<ngroups; j++){
                    fRlm[avg][j][i] = (double)REAL(poppars)[j*(nyears-1)+i];
                    fRlm[var][j][i] = (double)REAL(poppars)[ngroups*(nyears-1)+j*(nyears-1)+i];
                }
            }
            len0[avg][i] = (double)REAL(poppars)[2*ngroups*(nyears-1)+i];
            len0[var][i] = (double)REAL(poppars)[2*ngroups*(nyears-1) + nyears +i];
        }
        year0[avg][0] = (double)REAL(poppars)[2*ngroups*(nyears-1)];
        year0[var][0] = (double)REAL(poppars)[2*ngroups*(nyears-1) + nyears];
        for(int i=1;i<nages;i++){
            year0[avg][i] = (double)REAL(poppars)[2*ngroups*(nyears-1) + 2*nyears + i-1];
            year0[var][i] = (double)REAL(poppars)[2*ngroups*(nyears-1) + 2*nyears + nages-1 + i-1];
        }
        rB = (double)REAL(poppars)[2*ngroups*(nyears-1) + 2*nyears + 2*(nages-1) ];

    }
    else error("model indiciation incorrect (C code)\n");


    //make output file and pointers
    SEXP predicted_length = PROTECT(allocVector(VECSXP, 2));
    SET_VECTOR_ELT(predicted_length, avg, PROTECT(allocVector(VECSXP, nyears)) );
    SET_VECTOR_ELT(predicted_length, var, PROTECT(allocVector(VECSXP, nyears)) );

    double*** p_predicted_length = (double***) malloc(2* sizeof(double**));
    p_predicted_length[avg] = (double**) malloc(nyears* sizeof(double*));
    p_predicted_length[var] = (double**) malloc(nyears* sizeof(double*));

    for(int i=0; i<nyears;i++){
        SET_VECTOR_ELT(VECTOR_ELT(predicted_length, avg),i, PROTECT(allocVector(REALSXP, nages)));
        SET_VECTOR_ELT(VECTOR_ELT(predicted_length, var),i, PROTECT(allocVector(REALSXP, nages)));

        p_predicted_length[avg][i] = (double*) REAL(VECTOR_ELT(VECTOR_ELT(predicted_length, avg),i));
        p_predicted_length[var][i] = (double*) REAL(VECTOR_ELT(VECTOR_ELT(predicted_length, var),i));

    }

    // copy predictions from first year
    for(int j=0;j<nages;j++){
        p_predicted_length[avg][0][j] = year0[avg][j];
        p_predicted_length[var][0][j] = year0[var][j];
    }

    // calculate predictions
    for(int i=1; i<nyears; i++){
        p_predicted_length[avg][i][0] = len0[avg][i];
        p_predicted_length[var][i][0] = len0[var][i];

        for(int j=1; j<nages; j++){
            int k=0;
            if(grouppar == 1 ){
                while( k < ngroups-1 && j > (double)REAL(feedingbounds)[k] ) k++;
            }
            else if(grouppar == 2){
                while( k < ngroups-1 && p_predicted_length[avg][i-1][j-1] > (double)REAL(feedingbounds)[k] ) k++;
            }
            
            p_predicted_length[avg][i][j] = predict_length_avg(p_predicted_length[avg][i-1][j-1], fRlm[avg][k][i-1], rB, schrinking);
            p_predicted_length[var][i][j] = predict_length_var(p_predicted_length[var][i-1][j-1], fRlm[var][k][i-1], rB);

        }
    }

    // free memory
                 
    for(int i=0; i<ngroups; i++){
            free(fRlm[avg][i]);
            free(fRlm[var][i]);
        }

    free(fRlm[avg]);
    free(fRlm[var]);
    free(len0[avg]);
    free(len0[var]);
    free(year0[avg]);
    free(year0[var]);

    free(fRlm);
    free(len0);
    free(year0);

    // unprotect
    UNPROTECT(2*nyears+3);

    return predicted_length;


}



SEXP minLL_growthcurve(SEXP poppars, SEXP globpars, SEXP feedingbounds, SEXP samples, SEXP modeltype){

    int avg = 0;
    int var = 1;

    //split parameters
    int nyears = INTEGER(globpars)[1];
    int minyear = INTEGER(globpars)[2];
    int minage = INTEGER(globpars)[4];
    int nobs = INTEGER(globpars)[5];

    //calculate predictions
    SEXP predicted_length;


    PROTECT(predicted_length = predict_growthcurve(poppars, globpars, modeltype, feedingbounds));

    //pointer to predicted values
    double*** p_predicted_length = malloc(2*sizeof(double**));
    p_predicted_length[avg] = (double**) malloc(nyears*sizeof(double*));
    p_predicted_length[var] = (double**) malloc(nyears*sizeof(double*));

    for(int i=0; i<nyears; i++){
        p_predicted_length[avg][i] = (double*) REAL(VECTOR_ELT(VECTOR_ELT(predicted_length, avg),i));
        p_predicted_length[var][i] = (double*) REAL(VECTOR_ELT(VECTOR_ELT(predicted_length, var),i));
    }

    //pointers to data
    int* pyear = (int*) INTEGER(VECTOR_ELT(samples, 0));
    int* page = (int*) INTEGER(VECTOR_ELT(samples, 1));
    double* psize = (double*) REAL(VECTOR_ELT(samples, 2));
    double* pweight = (double*) REAL(VECTOR_ELT(samples, 3));

    // calculate sum of log likelyhood
    double LL = 0;

    for(int i=0; i<nobs; i++){
        LL += pweight[i] * log( calculate_prob(psize[i], p_predicted_length[avg][pyear[i]-minyear][page[i]-minage], p_predicted_length[var][pyear[i]-minyear][page[i]-minage]) );
    }

    UNPROTECT(1);

    return ScalarReal(-LL);

}









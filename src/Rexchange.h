/*
 * Rexchange.h
 *
 *  Created on: 03.05.2009
 *  Changed on: 10.12.2012
 *      Author: daniel (adapted by andrea)
 */

#ifndef REXCHANGE_H_
#define REXCHANGE_H_

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


// interface to R via C
void _myHyp2f1(double *res, double *a, double *b, double *c, double *z, int *n);
//void myHyperg2F1(double *res, double *a, double *b, double *c, double *z, int *n);

// register functions for R:
R_NativePrimitiveArgType myArgs[6] = {REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
R_CMethodDef cMethods[] = 
{
    { "_myHyp2f1", (DL_FUNC) &_myHyp2f1, 6, myArgs},
    //{ "myHyperg2F1", (DL_FUNC)&myHyperg2F1, 6, myArgs},
    { NULL, NULL, 0 }    
};

void
R_init_Repitools(DllInfo *info)
{
    R_registerRoutines(info,
		       cMethods, // .C
		       NULL, // .Call
		       NULL, // .Fortran
		       NULL); // .External
}




#endif /* REXCHANGE_H_ */

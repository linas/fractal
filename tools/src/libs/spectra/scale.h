
/*
 * FUNCTION:
 * This module computes a measure of the scalability of stochastic
 * processes.
 * 
 * HISTORY:
 * Linas Vepstas December 1992
 */

#ifdef AIX221
#ifndef _C_func
#define _C_func
#endif _C_func
#endif

#ifdef AIX315
#define ANSI_C
#endif

#include <math.h>

/* function prototypes */

#ifdef ANSI_C

extern double DifferentialVariance (double *arr, int narr, double scale);
extern double DifferentialRMS (double *arr, int narr, double scale);
extern double DifferentialCoVariance (double *xarr, double *yarr, int narr, double scale);

#else

extern double DifferentialVariance (); 
extern double DifferentialRMS (); 
extern double DifferentialCoVariance (); 

#endif

/* ====================== END OF FILE ========================= */

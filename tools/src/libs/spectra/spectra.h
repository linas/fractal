
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

extern void Fourrier (double *arr, int narr, 
                    double omega, 
                    double *ao, double *bo);
extern double SpectralDensity (double *arr, int narr, double omega);

extern void ConvolveWithGaussian (double *arr, int narr, 
                           double sigma,
                           double **garr, int *gnarr);
extern void ConvolveFourrier (double *arr, int narr, 
                    double omega, double sigma,
                    double *ao, double *bo);
extern double ConvolveSpectralDensity (double *arr, int narr, 
                                double omega, double sigma);
extern void SmoothFourrier (double *arr, int narr, 
                    double omega, 
                    double *ao, double *bo);
extern double SmoothSpectralDensity (double *arr, int narr, double omega);


#else

extern void Fourrier ();
extern double SpectralDensity (); 
extern void ConvolveWithGaussian ();
extern void ConvolveFourrier ();
extern double ConvolveSpectralDensity ();
extern void SmoothFourrier ();
extern double SmoothSpectralDensity ();

#endif

/* ====================== END OF FILE ========================= */

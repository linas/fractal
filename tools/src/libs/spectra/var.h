
/*
 * variance and covariance
 *
 * FUNCTION:
 * This module contains routines that compute the variance and
 * covarience of arrays of data
 *
 * HISTORY:
 * Linas Vepstas December 1993
 */

#ifndef _C_func
#define _C_func
#endif _C_func

#include <math.h>

/* ======================================================== */
/* compute the mean in a data array */

#ifdef ANSI_C
double Mean (double *arr, int narr);
double Variance (double *arr, int narr);
double StandardDeviation (double *arr, int narr);
double CentralMoment (double *arr, int narr, int moment);
double CoVariance (double *xarr,  double *yarr, int narr);
double Correlation (double *xarr,  double *yarr, int narr);
double AutoCoVariance (double *xarr, int narr, int delta);
double AutoCorrelation (double *xarr, int narr, int delta);
double Stationarity (double *xarr, int narr, int delta);

double GMean (double *arr, int narr, double center, double width);
double GVariance (double *arr, int narr, double center, double width);
double GStandardDeviation (double *arr, int narr, double center, double width);
double GCentralMoment (double *arr, int narr, int moment, double center, double width);
double GCoVariance (double *xarr,  double *yarr, int narr, double center, double width);
double GCorrelation (double *xarr,  double *yarr, int narr, double center, double width);
double GAutoCoVariance (double *xarr, int narr, int delta, double center, double width);
double GAutoCorrelation (double *xarr, int narr, int delta, double center, double width);

#else
double Mean ();
double Variance ();
double StandardDeviation ();
double CentralMoment ();
double CoVariance ();
double Correlation ();
double AutoCoVariance ();
double AutoCorrelation ();
double Stationarity ();

double GMean ();
double GVariance ();
double GStandardDeviation ();
double GCentralMoment ();
double GCoVariance ();
double GCorrelation ();
double GAutoCoVariance ();
double GAutoCorrelation ();

#endif

/* ================== END OF FILE ========================= */


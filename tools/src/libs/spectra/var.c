
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
#include "var.h"


/* ======================================================== */
/* compute the mean in a data array */

#ifdef ANSI_C
double Mean (double *arr, int narr)
#else
double Mean (arr, narr)
double *arr;
int narr;
#endif
{
   int i;
   double ns, xmean;

   /* The mean and varience of this differential will be computed */
   xmean = 0.0;
   for (i=0; i<narr; i++) {
      xmean += arr[i];
   }

   /* normalize mean and variance */
   ns = (double) (narr);
   xmean /= ns;

   return (xmean);
}

/* ======================================================== */
/* compute the variaence in a data array */

#ifdef ANSI_C
double Variance (double *arr, int narr)
#else
double Variance (arr, narr)
double *arr;
int narr;
#endif
{
   int i;
   double ns, xmean, xvar;

   /* The mean and varience of this differential will be computed */
   xmean = 0.0;
   xvar = 0.0;
   for (i=0; i<narr; i++) {
      xmean += arr[i];
      xvar += arr[i]*arr[i];
   }

   /* normalize mean and variance */
   ns = (double) (narr);
   xmean /= ns;
   xvar /= ns;

   /* compute mean square deviation */
   xvar -= xmean*xmean;

   return (xvar);
}

/* ======================================================== */
/* compute the standard deviation in a data array */

#ifdef ANSI_C
double StandardDeviation (double *arr, int narr)
#else
double StandardDeviation (arr, narr)
double *arr;
int narr;
#endif
{
   return ( sqrt(Variance (arr, narr)));
}

/* ======================================================== */
/* compute the nth central moment in a data array */
/* warning -- this is slow -- it uses the pow math library function */

#ifdef ANSI_C
double CentralMoment (double *arr, int narr, int moment)
#else
double CentralMoment (arr, narr, moment)
double *arr;
int narr, moment;
#endif
{
   int i;
   double ns, xmean, xvar;

   /* The mean and varience of this differential will be computed */
   xmean = Mean (arr, narr);
   xvar = 0.0;
   for (i=0; i<narr; i++) {
      xvar += pow ( (arr[i] - xmean), (double) moment);
   }

   /* normalize the moment */
   ns = (double) (narr);
   xvar /= ns;

   return (xvar);
}

/* ======================================================== */
/* compute the co-variance in a data array */

#ifdef ANSI_C
double CoVariance (double *xarr,  double *yarr, int narr)
#else
double CoVariance (xarr, yarr, narr)
double *xarr, *yarr;
int narr;
#endif
{
   int i;
   double ns, xmean, ymean, cov;

   /* The mean and varience of this differential will be computed */
   xmean = 0.0;
   ymean = 0.0;
   cov = 0.0;
   for (i=0; i<narr; i++) {
      xmean += xarr[i];
      ymean += yarr[i];
      cov += xarr[i]*yarr[i];
   }

   /* normalize mean and variance */
   ns = (double) (narr);
   xmean /= ns;
   ymean /= ns;
   cov /= ns;

   /* compute mean square deviation */
   cov -= xmean*ymean;

   return (cov);
}

/* ======================================================== */
/* compute the correleation in a data array */
/* The correlation is defined as:
 * r = CoVariance (xarr, yarr) / sqrt (Variance (xarr) * Varience (yarr))
 * Note that r runs between +1 and -1
 *
 * The code below is a tad more computationally efficient for getting
 * correlation than the above formula would lead to.
 */

#ifdef ANSI_C
double Correlation (double *xarr,  double *yarr, int narr)
#else
double Correlation (xarr, yarr, narr)
double *xarr, *yarr;
int narr;
#endif
{
   int i;
   double ns, xmean, ymean, cov;
   double xvar, yvar;

   /* The mean and varience of this differential will be computed */
   xmean = 0.0;
   ymean = 0.0;
   xvar = 0.0;
   yvar = 0.0;
   cov = 0.0;
   for (i=0; i<narr; i++) {
      xmean += xarr[i];
      ymean += yarr[i];
      xvar += xarr[i] * xarr[i];
      yvar += yarr[i] * yarr[i];
      cov += xarr[i]*yarr[i];
   }

   /* normalize mean and variance */
   ns = (double) (narr);
   xmean /= ns;
   ymean /= ns;
   xvar /= ns;
   yvar /= ns;
   cov /= ns;

   /* compute mean square deviation */
   xvar -= xmean*xmean;
   yvar -= ymean*ymean;
   cov -= xmean*ymean;

   cov /= sqrt (xvar*yvar);

   return (cov);
}

/* ======================================================== */
/* compute the auto-co-variance in a data array */

#ifdef ANSI_C
double AutoCoVariance (double *xarr, int narr, int delta)
#else
double AutoCoVariance (xarr, narr, delta)
double *xarr;
int narr;
int delta;
#endif
{
   if (delta < 0) delta = - delta;
   if (delta >= narr) return (0.0);
   return (CoVariance (xarr, &xarr[delta], narr-delta));
}

/* ======================================================== */
/* compute the auto-correlation in a data array */

#ifdef ANSI_C
double AutoCorrelation (double *xarr, int narr, int delta)
#else
double AutoCorrelation (xarr, narr, delta)
double *xarr;
int narr;
int delta;
#endif
{
   if (delta < 0) delta = - delta;
   if (delta >= narr) return (0.0);
   return (Correlation (xarr, &xarr[delta], narr-delta));
}

/* ======================================================== */
/* compute the variance of the stationary increments of a process */

#ifdef ANSI_C
double Stationarity (double *xarr, int narr, int delta)
#else
double Stationarity (xarr, narr, delta)
double *xarr;
int narr;
int delta;
#endif
{
   double station;
   double *diff;
   int i;

   if (delta < 0) delta = - delta;
   if (delta >= narr) return (0.0);
   
   diff = (double *) malloc ((narr-delta) * sizeof (double));
   for (i=0; i<narr-delta; i++) {
      diff[i] = xarr[i] - xarr [i+delta];
   }
   station = Variance (diff, narr-delta);
   free (diff);

   return (station);
}

/* ======================================================== */
/* ======================================================== */
/* compute the mean in a data array */
/* data is weighted by a gaussian */

#ifdef ANSI_C
double GMean (double *arr, int narr, double center, double width)
#else
double GMean (arr, narr, center, width)
double *arr;
int narr;
double center, width;
#endif
{
   int i;
   double ns, xmean;
   double fac;
   double weight;
   double norm;

   /* The mean and varience of this differential will be computed */
   fac = 1.0 / (2.0 * width * width);
   norm = 0.0;
   xmean = 0.0;
   for (i=0; i<narr; i++) {
      weight = exp (- fac * ((double) i - center) * ((double) i - center));
      norm += weight;
      xmean += arr[i] * weight;
   }

   /* normalize mean and variance */
   ns = 1.0 / norm;
   xmean *= ns;

   return (xmean);
}

/* ======================================================== */
/* compute the variaence in a data array */
/* data is weighted by a gaussian */

#ifdef ANSI_C
double GVariance (double *arr, int narr, double center, double width)
#else
double GVariance (arr, narr, center, width)
double *arr;
int narr;
double center, width;
#endif
{
   int i;
   double ns, xmean, xvar;
   double fac, weight;
   double norm;

   /* The mean and varience of this differential will be computed */
   fac = 1.0 / (2.0 * width * width);
   norm = 0.0;
   xmean = 0.0;
   xvar = 0.0;
   for (i=0; i<narr; i++) {
      weight = exp (- fac * ((double) i - center) * ((double) i - center));
      norm += weight;
      xmean += arr[i] * weight;
      xvar += arr[i]*arr[i] * weight;
   }

   /* normalize mean and variance */
   ns = 1.0 / norm;
   xmean *= ns;
   xvar *= ns;

   /* compute mean square deviation */
   xvar -= xmean*xmean;

   return (xvar);
}

/* ======================================================== */
/* compute the standard deviation in a data array */
/* data is weighted by a gaussian */

#ifdef ANSI_C
double GStandardDeviation (double *arr, int narr, double center, double width)
#else
double GStandardDeviation (arr, narr, center, width)
double *arr;
int narr;
double center, width;
#endif
{
   return ( sqrt(GVariance (arr, narr, center, width)));
}

/* ======================================================== */
/* compute the nth central moment in a data array */
/* warning -- this is slow -- it uses the pow math library function */
/* data is weighted by a gaussian */

#ifdef ANSI_C
double GCentralMoment (double *arr, int narr, int moment, double center, double width)
#else
double GCentralMoment (arr, narr, moment, center, width)
double *arr;
int narr, moment;
double center, width;
#endif
{
   int i;
   double ns, xmean, xvar;
   double fac, weight;
   double norm;

   /* The mean and varience of this differential will be computed */
   fac = 1.0 / (2.0 * width * width);
   norm = 0.0;
   xmean = GMean (arr, narr, center, width);
   xvar = 0.0;
   for (i=0; i<narr; i++) {
      weight = exp (- fac * ((double) i - center) * ((double) i - center));
      norm += weight;
      xvar += pow ( (arr[i] - xmean), (double) moment) * weight;
   }

   /* normalize the moment */
   ns = 1.0 / norm;
   xvar *= ns;

   return (xvar);
}

/* ======================================================== */
/* compute the co-variance in a data array */
/* data is weighted by a gaussian */

#ifdef ANSI_C
double GCoVariance (double *xarr,  double *yarr, int narr, double center, double width)
#else
double GCoVariance (xarr, yarr, narr, center, width)
double *xarr, *yarr;
int narr;
double center, width;
#endif
{
   int i;
   double ns, xmean, ymean, cov;
   double fac, weight;
   double norm;

   /* The mean and varience of this differential will be computed */
   fac = 1.0 / (2.0 * width * width);
   norm = 0.0;
   xmean = 0.0;
   ymean = 0.0;
   cov = 0.0;
   for (i=0; i<narr; i++) {
      weight = exp (- fac * ((double) i - center) * ((double) i - center));
      norm += weight;
      xmean += xarr[i] * weight;
      ymean += yarr[i] * weight;
      cov += xarr[i]*yarr[i] * weight;
   }

   /* normalize mean and variance */
   ns = 1.0 / norm;
   xmean *= ns;
   ymean *= ns;
   cov *= ns;

   /* compute mean square deviation */
   cov -= xmean*ymean;

   return (cov);
}

/* ======================================================== */
/* compute the correleation in a data array */
/* The correlation is defined as:
 * r = CoVariance (xarr, yarr) / sqrt (Variance (xarr) * Varience (yarr))
 * Note that r runs between +1 and -1
 *
 * The code below is a tad more computationally efficient for getting
 * correlation than the above formula would lead to.
 */
/* data is weighted by a gaussian */

#ifdef ANSI_C
double GCorrelation (double *xarr,  double *yarr, int narr, double center, double width)
#else
double GCorrelation (xarr, yarr, narr, center, width)
double *xarr, *yarr;
int narr;
double center, width;
#endif
{
   int i;
   double ns, xmean, ymean, cov;
   double xvar, yvar;
   double fac, weight;
   double norm;

   /* The mean and varience of this differential will be computed */
   fac = 1.0 / (2.0 * width * width);
   norm = 0.0;
   xmean = 0.0;
   ymean = 0.0;
   xvar = 0.0;
   yvar = 0.0;
   cov = 0.0;
   for (i=0; i<narr; i++) {
      weight = exp (- fac * ((double) i - center) * ((double) i - center));
      norm += weight;
      xmean += xarr[i] * weight;
      ymean += yarr[i] * weight;
      xvar += xarr[i] * xarr[i] * weight;
      yvar += yarr[i] * yarr[i] * weight;
      cov += xarr[i]*yarr[i] * weight;
   }

   /* normalize mean and variance */
   ns = 1.0 / norm;
   xmean *= ns;
   ymean *= ns;
   xvar *= ns;
   yvar *= ns;
   cov *= ns;

   /* compute mean square deviation */
   xvar -= xmean*xmean;
   yvar -= ymean*ymean;
   cov -= xmean*ymean;

   cov /= sqrt (xvar*yvar);

   return (cov);
}

/* ======================================================== */
/* compute the auto-co-variance in a data array */

#ifdef ANSI_C
double GAutoCoVariance (double *xarr, int narr, int delta, double center, double width)
#else
double GAutoCoVariance (xarr, narr, delta, center, width)
double *xarr;
int narr;
int delta;
double center, width;
#endif
{
   if (delta < 0) delta = - delta;
   if (delta >= narr) return (0.0);
   return (GCoVariance (xarr, &xarr[delta], narr-delta, center, width));
}

/* ======================================================== */
/* compute the auto-correlation in a data array */

#ifdef ANSI_C
double GAutoCorrelation (double *xarr, int narr, int delta, double center, double width)
#else
double GAutoCorrelation (xarr, narr, delta, center, width)
double *xarr;
int narr;
int delta;
double center, width;
#endif
{
   if (delta < 0) delta = - delta;
   if (delta >= narr) return (0.0);
   return (GCorrelation (xarr, &xarr[delta], narr-delta, center, width));
}

/* ================== END OF FILE ========================= */


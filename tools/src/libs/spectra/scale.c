
#ifndef _C_func
#define _C_func
#endif _C_func

#include <math.h>
#include "scale.h"

/* ======================================================== */
/*
 * This module computes a measure of the scalability of stochastic
 * processes.
 * 
 * The measure is basically something funky I've invented.
 * 
 * The input data is smoothed with a gaussian of width "sigma".
 * Taking the first derivative of the smoothed data, we then
 * compute the RMS variation of this derivative.
 * 
 * Algorithmically, we do this the other way around -- by employing
 * integration by parts, we essentially convolve the first derivative of
 * the gaussian with the data, and take the RMS of the convolution.
 * 
 * This routine returns one floating point value -- the RMS of the 
 * first derivative of the smoothed data, for the given sigma factor.
 *
 * Linas Vepstas December 1992
 */

#ifdef ANSI_C
double DifferentialCoVariance (double *xarr, double *yarr, int narr, double sigma) 
#else
double DifferentialCoVariance (xarr, yarr, narr, sigma) 
double *xarr, *yarr, sigma;
int narr;
#endif
{
   double *coeff;
   int ncoeff;
   double x, y;
   int i, j;
   double norm, frac;
   double xmean, ymean, cov;
   double ns;

   /* malloc an array to store differential coefficients */
   /* the number of coefficients here is about sufficient to give 
      us 1.0e-10 accuracy */
   ncoeff = (int) (7.1*sigma);
   coeff = (double *) malloc (ncoeff*sizeof (double));

   /* initialize the coefficient array */
   /* Compute a normalized gaussian */
   for (i=0; i<ncoeff; i++) {
      x = (double) i;
      y = - x*x / (sigma*sigma);
      y = exp (y);
      coeff[i] = y;
   }

   /* normalize */
   norm = 0.0;
   for (i=1; i<ncoeff; i++) {
      norm += coeff[i];
   }
   norm *= 2.0;
   norm += coeff[0];
   norm = 1.0 / norm;
   
   for (i=0; i<ncoeff; i++) coeff[i] *= norm;

   /* take finite differences */
   for (i=ncoeff-1; i>0; i--) coeff[i] -= coeff[i-1];
   coeff [0] = 0.0;

   /* correct round-off errors for large sigmas */
   ns = (int) (0.06 * sigma);
   for (i=1; i<ns; i++) {
      x = (double) i;
      y = - x*x / (sigma*sigma);
      /* taylor expansion of exp(y) to sixth power */
      coeff[i] = 1.0 + y/6.0;
      coeff[i] = 1.0 + 0.2 * y * coeff[i];
      coeff[i] = 1.0 + 0.25 * y * coeff[i];
      coeff[i] = 1.0 + y * coeff[i] / 3.0;
      coeff[i] = 1.0 + 0.5 * y * coeff[i];
      coeff[i] = norm * y * coeff[i];
   }
   for (i=ns-1; i>0; i--) coeff[i] -= coeff[i-1];

   /* The mean and varience of this differential will be computed */
   xmean = 0.0;
   ymean = 0.0;
   cov = 0.0;
   for (i=ncoeff-2; i<narr-ncoeff+1; i++) {

      /* convolve differential coefficients with data */
      x = 0.0;
      y = 0.0;
      for (j=1; j<ncoeff; j++) {
         x +=  (xarr[i-j+1] - xarr [i+j]) * coeff[j];
         y +=  (yarr[i-j+1] - yarr [i+j]) * coeff[j];
      }

      xmean += x;
      ymean += y;
      cov += x*y;
   }

   /* normalize mean and variance */
   ns = (double) (narr - 2*ncoeff + 3);
   xmean /= ns;
   ymean /= ns;
   cov /= ns;

   /* compute the variance */
   cov -= xmean*ymean;

   /* rescale rms deviation */
   cov *= sigma;

   /* free up allocated memory */
   free (coeff);

   return (cov);
}

/* ======================================================== */
/*
 * This module computes a measure of the scalability of stochastic
 * processes.
 * 
 * The measure is basically something funky I've invented.
 * 
 * The input data is smoothed with a gaussian of width "sigma".
 * Taking the first derivative of the smoothed data, we then
 * compute the RMS variation of this derivative.
 * 
 * Algorithmically, we do this the other way around -- by employing
 * integration by parts, we essentially convolve the first derivative of
 * the gaussian with the data, and take the RMS of the convolution.
 * Painful experiance shows, however, that rather than using the
 * ananlytic derivative of the gaussian, the finite difference should be
 * used instead. Turns out the anaylytic derivative is a bad approximation
 * when the gaussian is narrow (there are insufficient data points).
 * 
 * This routine returns one floating point value -- the RMS of the 
 * first derivative of the smoothed data, for the given sigma factor.
 *
 * Linas Vepstas December 1992
 */

#ifdef ANSI_C
double DifferentialVariance (double *arr, int narr, double sigma) 
#else
double DifferentialVariance (arr, narr, sigma) 
double *arr, sigma;
int narr;
#endif
{
   double *coeff;
   int ncoeff;
   double x, y;
   int i, j;
   double norm, frac;
   double xmean, xvar;
   double ns;

   /* malloc an array to store differential coefficients */
   /* the number of coefficients here is about sufficient to give 
      us 1.0e-10 accuracy */
   ncoeff = (int) (7.1*sigma);
   coeff = (double *) malloc (ncoeff*sizeof (double));

#define DERIV_GAUSS
#ifdef DERIV_GAUSS
   /* initialize the coefficient array */
   /* The coeffiecients are the first derivative of a normalized
    * gaussian. */

#ifdef ANALYT_DERIV
/*
Warning -- don't use this -- it give bad behaviour for narrow
gaussians.
*/
   for (i=0; i<ncoeff; i++) {
      x = (double) i;
      y = - x*x / (sigma*sigma);
      y = exp (y);
      coeff[i] = x * y;
   }

   /* The normalization cannot be gotten from an exact expression, but
    * must be computed. The expression used hare can be derived by
    * asking that the points on the gaussian sum to one. */
   /* a = 2.0 / (sigma*sigma*sigma * sqrt (M_PI)); */
   norm = 0.0;
   for (i=1; i<ncoeff; i++) {
      frac = 0.0;
      for (j=i+1; j<ncoeff; j++) frac += coeff[j];
      frac *= 2.0;
      norm += coeff[i] + frac;
   }
   norm = -1.0 / norm;
   
   for (i=0; i<ncoeff; i++) coeff[i] *= norm;

#else /* ANALYT_DERIV */

   /* initialize the coefficient array */
   /* Compute a normalized gaussian */
   for (i=0; i<ncoeff; i++) {
      x = (double) i;
      y = - x*x / (sigma*sigma);
      y = exp (y);
      coeff[i] = y;
   }

   /* normalize */
   norm = 0.0;
   for (i=1; i<ncoeff; i++) {
      norm += coeff[i];
   }
   norm *= 2.0;
   norm += coeff[0];
   norm = 1.0 / norm;
   
   for (i=0; i<ncoeff; i++) coeff[i] *= norm;

   /* take finite differences */
   for (i=ncoeff-1; i>0; i--) coeff[i] -= coeff[i-1];
   coeff [0] = 0.0;

   /* correct round-off errors for large sigmas */
   ns = (int) (0.06 * sigma);
   for (i=1; i<ns; i++) {
      x = (double) i;
      y = - x*x / (sigma*sigma);
      /* taylor expansion of exp(y) to sixth power */
      coeff[i] = 1.0 + y/6.0;
      coeff[i] = 1.0 + 0.2 * y * coeff[i];
      coeff[i] = 1.0 + 0.25 * y * coeff[i];
      coeff[i] = 1.0 + y * coeff[i] / 3.0;
      coeff[i] = 1.0 + 0.5 * y * coeff[i];
      coeff[i] = norm * y * coeff[i];
   }
   for (i=ns-1; i>0; i--) coeff[i] -= coeff[i-1];

#endif /* ANALYT_DERIV */

   /* The mean and varience of this differential will be computed */
   xmean = 0.0;
   xvar = 0.0;
   for (i=ncoeff-2; i<narr-ncoeff+1; i++) {

      /* convolve differential coefficients with data */
      x = 0.0;
      for (j=1; j<ncoeff; j++) {
         x +=  (arr[i-j+1] - arr [i+j]) * coeff[j];
      }
      xmean += x;
      xvar += x*x;
   }

   /* normalize mean and variance */
   ns = (double) (narr - 2*ncoeff + 3);

#endif /* DERIV_GAUSS */

#ifdef STRAIGHT_GAUSS

#ifdef INTEG_GAUSS
/*
Warning -- don't use this -- it give bad behaviour for narrow
gaussians.
*/
   /* initialize the coefficient array */
   /* The coeffiecients are the first derivative of a normalized
    * gaussian */
   frac = 0.0;
   for (i=0; i<ncoeff; i++) {
      x = (double) i;
      y = - x*x / (sigma*sigma);
      y = exp (y);
      coeff[i] = -x * y;
      frac += coeff[i];
   }

   /* integrate */
   coeff [0] = - frac;
   for (i=1; i<ncoeff; i++) {
      coeff[i] += coeff[i-1];
   }
#else /* INTEG_GAUSS */

   /* initialize the coefficient array */
   /* The coeffiecients are a normalized gaussian */
   for (i=0; i<ncoeff; i++) {
      x = (double) i;
      y = - x*x / (sigma*sigma);
      y = exp (y);
      coeff[i] = y;
   }

#endif /* INTEG_GAUSS */

   /* The normalization cannot be gotten from an exact expression, but
    * must be computed. The expression used hare can be derived by
    * asking that the points on the gaussian sum to one. */
   norm = 0.0;
   for (i=1; i<ncoeff; i++) {
      norm += coeff[i];
   }
   norm *= 2.0;
   norm += coeff[0];
   norm = 1.0 / norm;
   
   for (i=0; i<ncoeff; i++) coeff[i] *= norm;

   /* The mean and varience of this differential will be computed */
   xmean = 0.0;
   xvar = 0.0;
   for (i=ncoeff-1; i<narr-ncoeff; i++) {

      /* convolve differential coefficients with data */
      x =  (arr[i+1] - arr [i]) * coeff[0];
      for (j=1; j<ncoeff; j++) {
         x +=  (arr[i+j+1] - arr [i+j]) * coeff[j];
         x +=  (arr[i-j+1] - arr [i-j]) * coeff[j];
      }
      xmean += x;
      xvar += x*x;
   }

   /* normalize mean and variance */
   ns = (double) (narr - 2*ncoeff + 1);
#endif /* STRAIGHT_GAUSS */

   xmean /= ns;
   xvar /= ns;

   /* compute the variance */
   xvar -= xmean*xmean;

   /* rescale rms deviation */
   xvar *= sigma;

   /* free up allocated memory */
   free (coeff);

   return (xvar);
}

/* ======================================================== */

#ifdef JUNK
double Gag (arr, narr, sigma) 
double *arr, sigma;
int narr;
{
   double *coeff;
   int ncoeff;
   double x, y;
   int i, j;
   double norm, frac;
   double xmean, xvar;
   double ns;

   /* malloc an array to store differential coefficients */
   /* the number of coefficients here is about sufficient to give 
      us 1.0e-10 accuracy */
   ncoeff = (int) (7.1*sigma);
   coeff = (double *) malloc (ncoeff*sizeof (double));

   /* initialize the coefficient array */
   /* The coeffiecients are a normalized gaussian */
   for (i=0; i<ncoeff; i++) {
      x = (double) i;
      y = - x*x / (sigma*sigma);
      y = exp (y);
      coeff[i] = y;
   }

   /* The normalization cannot be gotten from an exact expression, but
    * must be computed. The expression used hare can be derived by
    * asking that the points on the gaussian sum to one. */
   norm = 0.0;
   for (i=1; i<ncoeff; i++) {
      norm += coeff[i];
   }
   norm *= 2.0;
   norm += coeff[0];
   norm = 1.0 / norm;
   
   for (i=0; i<ncoeff; i++) coeff[i] *= norm;

for (i=0; i<ncoeff; i++) printf ("its %d %g \n", coeff[i]);

   /* The mean and varience of this differential will be computed */
   xmean = 0.0;
   xvar = 0.0;
   for (i=ncoeff-1; i<narr-ncoeff; i++) {

      /* convolve differential coefficients with data */
      x =  (arr[i+1] - arr [i]) * coeff[0];
      for (j=1; j<ncoeff; j++) {
         x +=  (arr[i+j+1] - arr [i+j]) * coeff[j];
         x +=  (arr[i-j+1] - arr [i-j]) * coeff[j];
      }
{
double last;
j=ncoeff-1;
last = (arr[i+j+1] - arr [i+j]) * coeff[j];
last += (arr[i-j+1] - arr [i-j]) * coeff[j];
last /= x;
printf ("last erm is %g \n", last);
}
      xmean += x;
      xvar += x*x;
   }

   /* normalize mean and variance */
   ns = (double) (narr - 2*ncoeff +1);
   xmean /= ns;
   xvar /= ns;
printf ("straight %f %f ", xmean, xvar);
fflush (stdout);

   /* take finite differences */
   for (i=ncoeff-1; i>0; i--) coeff[i] -= coeff[i-1];
for (i=0; i<ncoeff; i++) printf ("its %d %g \n", coeff[i]);
   coeff [0] = 0.0;

   /* The mean and varience of this differential will be computed */
   xmean = 0.0;
   xvar = 0.0;
   for (i=ncoeff-2; i<narr-ncoeff+1; i++) {

      /* convolve differential coefficients with data */
      x = 0.0;
      for (j=1; j<ncoeff; j++) {
         x +=  (arr[i-j+1] - arr [i+j]) * coeff[j];
      }
{
double last;
j=ncoeff-1;
last = (arr[i-j+1] - arr [i+j]) * coeff[j];
last /= x;
printf ("last baerm is %g \n", last);
}
      xmean += x;
      xvar += x*x;
   }

   /* normalize mean and variance */
   ns = (double) (narr - 2*ncoeff +3);
   xmean /= ns;
   xvar /= ns;
printf ("finite %f %f \n", xmean, xvar);




   /* compute the variance */
   xvar -= xmean*xmean;

   /* rescale rms deviation */
   xvar *= sigma;

   /* free up allocated memory */
   free (coeff);

   return (xvar);
}
#endif /* JUNK */

/* ======================================================== */
/* rms differential variance */

#ifdef ANSI_C
double DifferentialRMS (double *arr, int narr, double sigma) 
#else
double DifferentialRMS (arr, narr, sigma) 
double *arr, sigma;
int narr;
#endif
{
   return (sqrt (DifferentialVariance(arr, narr, sigma)));
}

/* ================== END OF FILE ========================= */


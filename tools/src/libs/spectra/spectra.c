
#ifndef _C_func
#define _C_func
#endif _C_func

#include <math.h>
#include "spectra.h"

#include <stdio.h>

/* ======================================================== */
/*
 * Compute Fourrier transform  coefficients of data
 *
 * Input parameters:
 * double *arr: input array of real-valued data points
 * int narr:    number of elements in array
 * double omega:  input frequency. 
 * 
 * Output parameters:
 * double *ao:  sine fourier amplitude
 * double *bo:  cosine fourier amplitude
 *
 * The fourrier amplitudes are computed in the following manner:
 * ao = (1/N) sum (n=0 to N) arr[n] sine (omega*n)
 * bo = (1/N) sum (n=0 to N) arr[n] cosine (omega*n)
 * 
 */

#ifdef ANSI_C
void Fourrier (double *arr, int narr, double omega, double *ao, double *bo)
#else
void Fourrier (arr, narr, omega, ao, bo)
double *arr, omega, *ao, *bo;
int narr;
#endif
{
   double s, c, sn, cn, sm, cm;
   int i, j;
   double x, a, b;


   /* initialize sine and cosine coefficients */
   sn = 0.0;
   cn = 1.0;
   x = omega;
   s = sin (x);
   c = cos (x);

   /* integrate */
   a = 0.0;
   b = 0.0;
   for (i=0; i<narr; i++) {
      a += arr[i] * sn;
      b += arr[i] * cn;

      /* compute next sine and cosine values */
      sm = sn;
      cm = cn;
      sn = s*cm + c*sm;
      cn = c*cm - s*sm;
   }

   /* normalize */
   a /= ((double) narr);
   b /= ((double) narr);

   /* a is sine, b is cosine */
   *ao = a;
   *bo = b;

}

/* ======================================================== */
/*
 * Compute the spectral density of an array of data
 * The spectral density is, of course, just the sum of the squares of
 * the fourier components.
 */

#ifdef ANSI_C
double SpectralDensity (double *arr, int narr, double omega) 
#else
double SpectralDensity (arr, narr, omega)
double *arr, omega;
int narr;
#endif
{
   double x, si, co;

   Fourrier (arr, narr, omega, &si, &co);

   x = si*si + co*co;

   return (x);

}

/* ======================================================== */
/*
 * Compute Convolution of data with Gaussian
 *
 * Input parameters:
 * double *arr: input array of real-valued data points
 * int narr:    number of elements in array
 * double sigma:  width of gaussian
 * 
 * Output parameters:
 * double **garr:  array containing convolved data
 * int *gnarr:  size of convoluted data array
 *
 * The convolution is done in the following manner:
 * garr[i] = N * sum (n=0 to N-1) arr[n] * exp ( -(i-n)^2 / sigma^2)
 *
 * where N is a normalization factor:
 * N = sum (i = -INF to +INF) exp ( - i^2 / sigma^2)
 *
 */

#ifdef ANSI_C
void ConvolveWithGaussian (double *arr, int narr, 
                           double sigma,
                           double **garr, int *gnarr)
#else
void ConvolveWithGaussian (arr, narr, sigma, garr, gnarr)
double *arr, sigma, **garr;
int narr, *gnarr;
#endif
{
   double *coeff;
   int ncoeff;
   int i, j, ns;
   double x, y, norm;

   /* malloc a coefficient array */
   ncoeff = (int) (5.1 * sigma) + 1;
   *gnarr = narr - 2*ncoeff+2;
   
   /* if not enough data, don't smooth */
   if ((*gnarr <=0) || (ncoeff <2)) {
      *garr = (double *) malloc (narr * sizeof (double));
      for (i=0; i<narr; i++) (*garr)[i] = arr[i];
      *gnarr = narr;
      return;
   }

   coeff = (double *) malloc (ncoeff * sizeof (double));
   if (coeff == (double *) 0x0) {
      fprintf (stderr, "ConvolveWithGaussian: Unable to malloc \n");
      exit(1);
   }

   /* compute gaussian coefficients */
   for (i=0; i<ncoeff; i++) {
      x = (double) i;
      y = - x*x / (sigma*sigma);
      y = exp (y);
      coeff[i] = y;
   }

   /* normalize the gaussian */
   norm = 0.0;
   for (i=1; i<ncoeff; i++) {
      norm += coeff[i];
   }
   norm *= 2.0;
   norm += coeff[0];
   norm = 1.0 / norm;
   for (i=0; i<ncoeff; i++) coeff[i] *= norm;

   *garr = (double *) malloc ((*gnarr) * sizeof(double));
   if ((*garr) == (double *) 0x0) {
      fprintf (stderr, "ConvolveWithGaussian: Unable to malloc \n");
      exit(1);
   }

   /* convolve the raw data with the gaussian */
   for (i=ncoeff-1; i<narr-ncoeff+1; i++) {

      /* convolve differential coefficients with data */
      x = arr[i] * coeff[0];
      for (j=1; j<ncoeff; j++) {

         x +=  (arr[i+j] + arr [i-j]) * coeff[j];
      }
      (*garr)[i-ncoeff+1] = x;
   }

   free (coeff);
   return;
}

/* ======================================================== */
/*
 * Compute Smoothed Fourrier transform  coefficients of data
 *
 * Input parameters:
 * double *arr: input array of real-valued data points
 * int narr:    number of elements in array
 * double omega:  input frequency. 
 * double sigma:  width of gaussian
 * 
 * Output parameters:
 * double *ao:  sine fourier amplitude
 * double *bo:  cosine fourier amplitude
 *
 * The fourrier amplitudes are computed in the following manner:
 * ao = (1/N) sum (n=0 to N) arr[n] sine (omega*n)
 * bo = (1/N) sum (n=0 to N) arr[n] cosine (omega*n)
 *
 * The data is first smoothed with a gaussian. The gaussian is choosen
 * to be very narrow, to that it's effect is less than 1 part in 1e10.
 * (i.e. the smoothing wipes out the high frequencies, but changes the
 * fundamental by less than 1e-10).
 * 
 */

#ifdef ANSI_C
void ConvolveFourrier (double *arr, int narr, 
                    double omega, double sigma,
                    double *ao, double *bo)
#else
void ConvolveFourrier (arr, narr, omega, sigma, ao, bo)
double *arr, omega, sigma, *ao, *bo;
int narr;
#endif
{
   double norm;
   double *sarr;
   int ns;

   ConvolveWithGaussian (arr, narr, sigma, &sarr, &ns);

   /* compute the  Fourrier of the smoothed data */
   Fourrier (sarr, ns, omega, ao, bo);

   /* renormalize for convolution with gaussian */
   norm = 0.5 * sigma * omega;
   norm = exp (norm*norm);
   *ao *= norm;
   *bo *= norm;

   free (sarr);
   return;
}

/* ======================================================== */
/*
 * Compute the spectral density of an array of data
 * The spectral density is, of course, just the sum of the squares of
 * the fourier components.
 */

#ifdef ANSI_C
double ConvolveSpectralDensity (double *arr, int narr, double omega, double sigma) 
#else
double ConvolveSpectralDensity (arr, narr, omega, sigma)
double *arr, omega, sigma;
int narr;
#endif
{
   double x, si, co;

printf (" yo data sent at 0x%x \n", arr);
   ConvolveFourrier (arr, narr, omega, sigma, &si, &co);
   x = si*si + co*co;
   return (x);
}

/* ======================================================== */
/*
 * Compute Smoothed Fourrier transform  coefficients of data
 *
 * Input parameters:
 * double *arr: input array of real-valued data points
 * int narr:    number of elements in array
 * double omega:  input frequency. 
 * double sigma:  width of gaussian
 * 
 * Output parameters:
 * double *ao:  sine fourier amplitude
 * double *bo:  cosine fourier amplitude
 *
 * The fourrier amplitudes are computed in the following manner:
 * ao = (1/N) sum (n=0 to N) arr[n] sine (omega*n)
 * bo = (1/N) sum (n=0 to N) arr[n] cosine (omega*n)
 *
 * The data is first smoothed with a gaussian. The gaussian is choosen
 * to be very narrow, to that it's effect is less than 1 part in 1e10.
 * (i.e. the smoothing wipes out the high frequencies, but changes the
 * fundamental by less than 1e-10).
 * 
 */

#ifdef ANSI_C
void SmoothFourrier (double *arr, int narr, 
                    double omega, 
                    double *ao, double *bo)
#else
void SmoothFourrier (arr, narr, omega, ao, bo)
double *arr, omega, *ao, *bo;
int narr;
#endif
{
   ConvolveFourrier (arr, narr, omega, (1.0/omega), ao, bo);
   return;
}

/* ======================================================== */
/*
 * Compute the spectral density of an array of data
 * The spectral density is, of course, just the sum of the squares of
 * the fourier components.
 */

#ifdef ANSI_C
double SmoothSpectralDensity (double *arr, int narr, double omega) 
#else
double SmoothSpectralDensity (arr, narr, omega)
double *arr, omega;
int narr;
#endif
{
   double x, si, co;

   SmoothFourrier (arr, narr, omega, &si, &co);
   x = si*si + co*co;
   return (x);
}

/* ================== END OF FILE ========================= */


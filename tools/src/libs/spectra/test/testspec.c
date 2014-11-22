
#include <stdio.h>
#include "stoch.h"
#include "spectra.h"
#include <math.h>

/*
 * This module computes a measure of the scalability of stochastic
 * processes.
 * 
 * The measure is basically something funky I've invented.
 * Essentially, I compute the RMS variation of the first "derivative" of
 * the stochastic process.  The derivative is, in fact, not the point
 * derivative, but a smeared derivative, with the scale factor representing
 * the amount of smearing.  (The larger the scale factor, the more smeared
 * the derivative).
 *
 * Linas Vepstas December 1992
 */



#define NPOINTS 1000
main ()
{
   double data[NPOINTS];
   int i;
   double x, omega;

   StochInit ();

#ifdef WHITE
   /* generate a gaussian process */
   for (i=0; i<NPOINTS; i++) {
      data [i] = WhiteNoise ();
/*
      printf ("%d %f \n", i, data[i]);
*/
   }
#endif

#ifdef GAUSS
   /* generate a gaussian process */
   for (i=0; i<NPOINTS; i++) {
      data [i] = GaussianNoise ();
/*
      printf ("%d %f \n", i, data[i]);
*/
   }
#endif

#define BROWN
#ifdef BROWN
   /* generate a gaussian process */
   data[0] = GaussianNoise ();
   for (i=1; i<NPOINTS; i++) {
      data [i] = data[i-1] + GaussianNoise ();
/*
      printf ("%d %f \n", i, data[i]);
*/
   }
#endif

#define SPEC
#ifdef SPEC

#define NSAMP 1000
   for (i=1; i<=NSAMP; i++) {
/*
      x = DifferentialAutoVariance (data, NPOINTS, (double) i);
*/
      omega =((double) i) * M_PI / ((double) NSAMP);
      x = SpectralDensity (data, NPOINTS, omega);
/*
      x *= omega*omega;
*/
      x = log10 (x);
      omega = log10 (omega);

      printf ("%f %f \n", omega, x);
      fflush (stdout);
   }
#endif

}

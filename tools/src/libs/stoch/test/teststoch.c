
/*
 * FUNCTION:
 * simple test case for random number generator
 *
 * HISTORY:
 * Lians Vepstas December 1992
 */

#include <math.h>
#include "stoch.h"

main (argc, argv) 
int argc;
char * argv[];
{
   int i, n, num_samp;
   double x;
   double xm, xsq;

   if (argc != 2) {
      printf ("Usage: %s <num_samples> \n", argv[0]);
      exit (1);
   }
   num_samp = atoi (argv[1]);

   StochInit ();

   xm = xsq = 0.0;
   for (i=0; i< num_samp; i++) {
      x = GaussianNoise ();
/*
printf ("i %d g %g \n", i, xm);
*/
      xm += x;
      xsq += x*x;
   }
   
   xm /= (double) num_samp;
   xsq / = (double) num_samp;
   xsq -= xm*xm;
   xsq = sqrt (xsq);
   
   printf ("samples %d mean %g deviation %g invsqrt(nsamp) %g \n", 
      num_samp, xm, xsq, (1.0/sqrt ((double) num_samp)) );
  
   exit (0);
}

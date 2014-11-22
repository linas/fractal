
#include <sys/types.h>
#include <math.h>
#include "stoch.h"

main () {
   double sigma, H;
   double ref;

   int i;
   int nDeltas;
   double *delta;
   double fac, base;
   double x;

   StochInit ();

   sigma = 1.0;
   H = 0.5;

   /* H is the self-similarity parameter */
   /* in order to get approx 48 bits of precision, we need to recur down
    * this many levels */
   x = 48.0 / H;
   nDeltas = ((int) x) + 3;
   delta = (double *) malloc ( nDeltas * sizeof (double));

   fac = 1.0 - pow (2.0, (2.0*H-2.0));
   fac *= 0.5;
   fac = sigma * sqrt (fac);
   base = pow (0.5, H);

   /* initialize scaling factors */
   for (i=0; i<nDeltas; i++) {
      delta[i] = fac * pow (base, ((double) i));
      printf (" deltas are %d %g \n", i, delta[i]);
   }
   printf (" \n\n\n");

   ref = delta[0] * GaussianNoise ();
   GaussianNoise ();
   for (i=1; i<nDeltas; i++) {
      ref += delta[i] * GaussianNoise();
      GaussianNoise ();
      printf ("ref is %d %g \n", i, ref);
   }
   printf (" \n\n\n");

}

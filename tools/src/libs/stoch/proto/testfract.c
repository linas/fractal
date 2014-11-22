
/*
 * FUNCTION:
 * This program generates a fractal noise process
 *
 * HACK ALERT -- this program has some bug in it. the noise produced
 * contains discontinuites. OOPs.
 *
 * HISTORY:
 * created by Linas Vepstas December 1992
 */
#include <sys/types.h>
#include <stdio.h>
#include <math.h>
#include "stoch.h"

#ifdef DEBUG
#define Tpr(x) printf x;
#else
#define Tpr(x)
#endif

main (argc, argv)
int argc;
char **argv;
{
   double sigma, H;
   int n=0;

   int i;
   int nDeltas;
   double renorm;
   double *delta;
   double *sideA, *sideB;
   char *tree;
   double fac, base;
   double x;
   int nDepth;

   StochInit ();

   sigma = 1.0;
   if (argc < 2) {
      fprintf (stderr, "Usage: %s <H> \n", argv[0]);
      exit (1);
   }
   H = atof (argv[1]);
   

   /* H is the self-similarity parameter */
   /* in order to get approx 48 bits of precision, we need to recur down
    * this many levels */
   x = 40.0 / H;
   x = 10.0 / H;
   nDeltas = ((int) x);
   delta = (double *) malloc ( nDeltas * sizeof (double));
   sideA = (double *) malloc ( nDeltas * sizeof (double));
   sideB = (double *) malloc ( nDeltas * sizeof (double));
   tree = (char *) malloc ( nDeltas * sizeof (char));

   fac = 1.0 - pow (2.0, (2.0*H-2.0));
   fac *= 0.5;
   fac = sigma * sqrt (fac);
   base = pow (0.5, H);
   fac /= pow (base, ((double) nDeltas));

   /* initialize scaling factors */
   for (i=0; i<nDeltas; i++) {
      delta[i] = fac * pow (base, ((double) i));
      tree[i] = FALSE;
      Tpr ((" deltas are %d %g \n", i, delta[i]));
   }
   Tpr ((" \n\n\n"));

   /* set up midpoint recursion values */
   sideA [0] = delta[0] * GaussianNoise ();
   sideB [0] = delta[0] * GaussianNoise ();
   for (i=1; i<nDeltas; i++) {
      sideB [i] = (0.5 * (sideB[i-1] + sideA [i-1]));
      sideA [i] = sideA [i-1];
      sideA [i] += delta[i] * GaussianNoise();
      sideB [i] += delta[i] * GaussianNoise();
      Tpr (("sides are %d %g %g \n", i, sideA[i], sideB[i]));
   }
   Tpr ((" \n\n\n"));

   /* renormalize */
   renorm = sideA [nDeltas-1];
/*
   for (i=0; i<nDeltas; i++) {
      sideA [i] -= renorm;
      sideB [i] -= renorm;
   }
*/

   printf ("# \n");
   printf ("# Raw fractal noise data file \n");
   printf ("# h is the scaling parameter \n");
   printf ("h %f \n", H);
   
for (n=0; n<3000; n++) {

   /* perform binary addition, using carry bit */
   nDepth = nDeltas-1;
   for (i=nDeltas-1; i>=0; i--) {
      if (tree[i] > 0) {
         tree[i] = 0;
      } else {
         tree[i] = 1;
         nDepth = i;
         break;
      }
   }

   /* now, walk the binary tree, performing midpoint recursion */
   for (i=nDepth; i<nDeltas; i++) {
      if (tree[i]) {
         sideA [i] = sideB [i];
         sideB [i] = sideB [i-1];
         sideB [i] += delta[i] * GaussianNoise();
      } else {
         sideA [i] = sideA [i-1];
         sideB [i] = (0.5 * (sideB[i-1] + sideA[i-1]));
         sideA [i] += delta[i] * GaussianNoise();
         sideB [i] += delta[i] * GaussianNoise();
      }
      Tpr (("sides are %d %g %g \n", i, sideA[i], sideB[i]));
   }

/*
   printf (" %d %g \n", n, sideA[nDeltas-1]/renorm);
*/
   printf (" i %d f %g \n", n, sideA[nDeltas-1]);


   

}
   exit (0);
}

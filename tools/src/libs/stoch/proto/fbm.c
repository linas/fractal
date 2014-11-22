
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
   int npts, maxLevel;
   double renorm;
   double *delta, *pts;
   double fac, base;
   int d, D, N;

   int level;

   StochInit ();

   sigma = 1.0;

   if (argc < 2) {
      fprintf (stderr, "Usage: %s <H> \n", argv[0]);
      exit (1);
   }
   H = atof (argv[1]);
   
   maxLevel = 14;
   npts = 1;
   for (i=0; i<maxLevel; i++) npts *= 2;
   npts +=1;

   /* H is the self-similarity parameter */
   delta = (double *) malloc ( maxLevel * sizeof (double));
   pts = (double *) malloc ( npts * sizeof (double));

   fac = 1.0 - pow (2.0, (2.0*H-2.0));
   fac *= 0.5;
   fac = sigma * sqrt (fac);
   base = pow (0.5, H);

   /* initialize scaling factors */
   for (i=0; i<maxLevel; i++) {
      delta[i] = fac * pow (base, ((double) i));
      Tpr ((" deltas are %d %g \n", i, delta[i]));
   }
   Tpr ((" \n\n\n"));

   /* set up midpoint recursion values */
   N = npts-1;
   D = N;
   d = D/2;
   pts [0] = 0.0;
   pts [D] = delta[0]*GaussianNoise ();
   level = 1;
   while (level <= maxLevel) {
      for (i=d; i<N-d+1; i+=D) {
         pts[i] = 0.5 * (pts[i-d] + pts[i+d]);
      }
      for (i=0; i<N-d+1; i+=D) {
         pts[i] += delta[level] * GaussianNoise ();
      }
      D /= 2;
      d /= 2;
      level ++;
   }

   printf ("# \n");
   printf ("# Raw fractal noise data file \n");
   printf ("# h is the scaling parameter \n");
   printf ("h %f \n", H);

   for (i=0; i<npts; i++) {
      printf (" i %d f %g \n", i, pts[i]);
   }

   exit (0);
}

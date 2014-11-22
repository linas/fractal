
#include <sys/types.h>
#include <stdio.h>
#include <math.h>
#include "stoch.h"

#ifdef DEBUG
#define Tpr(x) printf x;
#else
#define Tpr(x)
#endif

/* This is algorithm InterpolatedFM1D from Peitgen-Saupe, page 88 */

void lacune (X, N, r, sigma, H)
double X[];
int N;
double r, sigma, H;
{
   double *Y;
   int mt, mT;
   double t, T, h;
   double or, delta;
   int i, idx;

   StochInit ();

   Y = (double *) malloc ((N+1) * sizeof (double));

   X[0] = 0.0;
   X[1] = sigma * GaussianNoise ();
   mT = 2;
   T = 1.0;
   or = 1.0 / r;

   while (mT < N) {
      mt = (int) (((double) mT) * or);
      if (mt == mT) mt= mT + 1;
      if (mt > N) mt = N;
   
      t = 1.0 / ((double) (mt-1));

      Y[0] = X[0];
      Y[mt-1] = X[mT-1];
      

      for (i=1; i<mt-1; i++) {
         idx = (int) ( ((double) i) * t / T);
         h =  ( ((double) i) * t / T) - ((double) idx);
         Y[i] = (1.0-h) * X[idx] + h * X[idx+1];
      }

      delta = sqrt (0.5) * pow (t, H) * sqrt (1.0 - pow (t/T, 2.0*(1.0-H)));

      for (i=0; i<mt; i++) {
         X[i] = Y[i] + delta * GaussianNoise ();
      }

      mT = mt;
      T = 1.0 / ((double) mT);
   }

   free (Y);
}


main (argc, argv)
int argc;
char **argv;
{
   double sigma, H;
   double r;
   int npts;
   int i;
   double *pts;

   int level;

   StochInit ();

   sigma = 1.0;
   r = 0.9;
   npts = 3000;

   if (argc < 2) {
      fprintf (stderr, "Usage: %s <H> \n", argv[0]);
      exit (1);
   }
   H = atof (argv[1]);

   pts = (double *) malloc ( (npts+1) * sizeof (double));

   lacune (pts, npts, r, sigma, H);

   for (i=0; i<npts; i++) {
      printf (" i %d f %g \n", i, pts[i]);
   }

}

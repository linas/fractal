
/*
 * walk.C
 *
 * random walk in theta
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

main () 
{
#define NBIN 100

   int i;
   int bin[NBIN];

   double gb0 = 1.0e-2;

   for (i=0; i<NBIN; i++)
   {
      bin[i] = 0;
   }

   double theta = 0.0;

   for (i=0; i<1230000; i++)
   {
      double x = rand() / ((double) RAND_MAX);
      double y = rand() / ((double) RAND_MAX);

      double d_theta = gb0 * x * sin(2.0*M_PI*y);
      theta += d_theta;
      if (0.0 > theta) { theta = -theta; }
      if (M_PI < theta) { theta = 2.0*M_PI-theta; }

      int ibin = (int) (((double)NBIN)*theta / M_PI);

      bin[ibin] ++;
   }

   for (i=0; i<NBIN; i++)
   {
      theta = (i+0.5)*M_PI / ((double) NBIN);
      printf ("%d	%g	%d\n", i, theta, bin[i]);
   }

}


/* 
 * FUNCTION:
 * Test operation of farey number converter
 *
 * HISTORY:
 * Linas Vepstas January 16 1994
 */

#include "Farey.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ------------------------------------------------------------ */

main (int argc, char *argv[])
{
   ContinuedFraction f;
   double xe, xo, t, g;
   int n, d;

   if (argc <3) {
      printf ("Usage: %s <num> <denom> \n", argv[0]);
      exit (1);
   }

   n = atoi (argv[1]);
   d = atoi (argv[2]);

   t = ((double) n) / ((double) d);
   /* f.SetReal (t); */
   f.SetRatio (n, d);

   xe = f.ToGapEven (); 
   xo = f.ToGapOdd (); 
	g = xe-xo;

	double w = 0.01;
	double xgap = f.ToXEven(w) - f.ToXOdd(w);

	g -= 1.0;
	double guess = (2*w - w*w + 0.5*w*w*w + 0.5*g*w*w*w)/((double)(d*d));

	double gg = xgap *d*d;
   gg /= w;
   gg -= 2.0;
	gg += w;
	gg -= 0.5* w*w;
	gg /= 0.5*w*w;

   printf ("\n------ gap(%d/%d) = %10.8g= %f - %f  gg=%g nterm=%d xg = %20.15g gus=%20.15g\n", 
						 n,d,g, xe, xo, gg, f.GetNumTerms(), xgap, guess);

}

/* ---------------------- END OF FILE ------------------------- */


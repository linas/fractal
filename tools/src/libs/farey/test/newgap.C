
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

	double guess = (2*w - w*w + 0.5*g*w*w*w)/((double)(d*d));

	g -= 1.0;
   printf ("\n------ gap(%d/%d) = %10.8g= %f - %f  nt=%d xg = %20.15g gus=%20.15g\n", 
						 n,d,g, xe, xo, f.GetNumTerms(), xgap, (xgap-guess)/xgap);

}

/* ---------------------- END OF FILE ------------------------- */


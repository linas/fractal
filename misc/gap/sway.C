
/*
 * sway.C
 *
 * simple graph of swaying
 *
 * Linas Vepstas October 2004
 */

#include "Farey.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


main (int argc, char *argv[])
{
	ContinuedFraction f;
	double w, t;
	int n,d;

   if (argc <2) {
      printf ("Usage: %s <denom> \n", argv[0]);
      exit (1);
   }

	d = atoi (argv[1]);
	printf ("#\n# denom=%d\n#\n", d);
  
	for (n=1; n<d; n++)
	{
		t = ((double) n) / ((double) d);
   	f.SetRatio (n,d);

		w = 1.0;
		double xm = f.ToXOdd(w);
		double xp = f.ToXEven(w);
		printf ("%g	%g	%g\n", t, xm, xp);
	}
}

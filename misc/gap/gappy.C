
/*
 * gappy.C
 *
 * Games with Gaps.
 *
 * Linas Vepstas October 2004
 */

#include "Farey.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NV 6

main (int argc, char *argv[])
{
	ContinuedFraction f;
	double w, t, gap[NV];
	int n,d,iw;

   if (argc <2) {
      printf ("Usage: %s <denom> \n", argv[0]);
      exit (1);
   }

	d = atoi (argv[1]);
	printf ("#\n# denom=%d\n#\n", d);
  
	for (n=0; n<=d; n++)
	{
		t = ((double) n) / ((double) d);
   	f.SetRatio (n,d);

		for (iw=0; iw<NV; iw++)
		{
			w = (double) (iw+1);
			w /= (double) NV;
			double xm = f.ToXOdd(w);
			double xp = f.ToXEven(w);
			gap[iw] = xp - xm;
		}
		printf ("%g	%g	%g	%g	%g	%g	%g\n", 
			t, gap[0], gap[1], gap[2], gap[3], gap[4], gap[5]);
	}
}

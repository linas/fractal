
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

/* return the greatest common factor */
int
gcf32 (int nume, int denom)
{
	int t;
	t = nume % denom;
	nume = denom;
	denom = t;

	/* Euclids algorithm for obtaining the gcf */
	while (0 != denom)
	{
		t = nume % denom;
		nume = denom;
		denom = t;
	}

	/* num now holds the GCD (Greatest Common Divisor) */
	return nume;
}

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
		int gcf, rn, rd;

		gcf = gcf32 (n,d);
		rn = n/gcf;
		rd = d/gcf;

		t = ((double) n) / ((double) d);
   	f.SetRatio (n,d);

		for (iw=0; iw<NV; iw++)
		{
			w = (double) (iw+1);
			w /= (double) NV;
			double xm = f.ToXOdd(w);
			double xp = f.ToXEven(w);
			gap[iw] = xp - xm;
			gap[iw] *= (double) rd*rd;
		}
		printf ("%g	%g	%g	%g	%g	%g	%g\n", 
			t, gap[0], gap[1], gap[2], gap[3], gap[4], gap[5]);
	}
}

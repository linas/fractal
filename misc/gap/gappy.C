
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

main (int argc, char *argv[])
{
	ContinuedFraction f;
	double w, t, xp, xm;
	int n,d;

	n = atoi (argv[1]);
	d = atoi (argv[2]);
  
	t = ((double) n) / ((double) d);
   f.SetRatio (n,d);

	for (w=0.0; w<2.0; w++)
	{
		xm = f.ToXMinus(w);
		xp = f.ToXPlus(w);
		printf ("its t=%g w=%g %g %g \n", t, w, xm, xp);
	}
}

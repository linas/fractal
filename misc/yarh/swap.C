
/* swap.C
 *
 * Swap terms in contineud fraction
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>

#include "Farey.h"

main (int argc, char *argv[])
{
	int i;
	ContinuedFraction f;
	// f.SetEvenize();
	
	if (argc<3)
	{
		printf ("Usage: %s <p> <q>\n", argv[0]);
		exit (1);
	}
	int p = atoi (argv[1]);
	int q = atoi (argv[2]);
	printf ("# swapping terms %d and %d\n", p,q);

	int nmax = 431;
	for (i=0; i<nmax; i++)
	{
		double x = ((double) i)/ ((double) (nmax));

		f.SetRatio (i,nmax);
		f.SwapTerms (p,q);
		double y = f.ToReal();

		printf ("%5d	%8.6g	%8.6g\n", i, x, y);

	}
}


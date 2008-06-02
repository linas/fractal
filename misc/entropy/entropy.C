/* 
 * entropy.C
 *
 * Explore the entropy of individual expansions of rational numbers.
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"

main (int argc, char *argv[])
{
	ContinuedFraction f;
	f.SetEvenize();

#if 0
	if (2>argc)
	{
		printf ("Usage: %s <z-value>\n", argv[0]);
		exit (1);
	}

	zre = atof (argv[1]);
#endif
	
	int nmax = 5531;
	int i;
	for (i=0; i<nmax; i++)
	{

		int p = i;
		int q = nmax;
		double x = ((double) p)/ ((double) q);
	
		f.SetRatio (p, q);

		int k;
		for (k=1; k<=f.GetNumTerms(); k++)
		{
			int j = f.GetTerm(k);
printf("duude k=%d j=%d\n", k, j);
		}
		
		double y;
		printf("%5d	%8.6g	%8.6g\n", i,x,y);

	}
}


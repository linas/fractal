/* 
 * defacto.C
 *
 * Explore the defacto entropy of rational numbers.
 *
 * Linas Vepstas October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gcf.h"
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

	int qmax;
	for (qmax = 2; qmax < 129; qmax++)
	{
		int total = 0;
		int totlen = 0;
	
		int p, q;
		for (q=1; q <= qmax; q++)
		{
			for (p=0; p<q; p++)
			{
				if (1 != gcf32(p,q)) continue;	
				f.SetRatio (p, q);

				int nterms = f.GetNumTerms();
				totlen += nterms;
				total ++;
			}
		}
	
		double avglen = ((double) totlen) / ((double) total);
		double ent = log(total) / M_LN2;
		ent /= avglen;
		// printf("%a total=%d avglen=%g ent=%g\n", total, avglen, ent);
		printf("%d	%d	%g	%g\n", qmax, total, avglen, ent);
	}
}


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

	unsigned long long total = 0LL;
	unsigned long long totlen = 0LL;
	
	int p, q;
	for (q=1; q <= 329000; q++)
	{
		for (p=0; p<q; p++)
		{
			if (1 != gcf32(p,q)) continue;	
			f.SetRatio (p, q);

			int nterms = f.GetNumTerms();
			totlen += nterms;
			total ++;
		}
	
		if ((q%100 == 0) || (q < 1000))
		{
			double avglen = ((double) totlen) / ((double) total);
			double ent = log((double)total) / M_LN2;
			ent /= avglen;
			// printf("%a total=%d avglen=%g ent=%g\n", total, avglen, ent);
			printf("%d	%lld	%g	%g\n", q, total, avglen, ent);
			fflush (stdout);
		}
	}
}


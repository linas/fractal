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

		int term[100];
		int cnt[100];
		int maxterm = 0;
		int nterms = f.GetNumTerms();
		int k;
		for (k=1; k<=nterms; k++)
		{
			int j = f.GetTerm(k);
			int notfound = 1;
			int p;
			for (p=0; p<maxterm; p++)
			{
				if (j == term[p]) {cnt[p] ++; notfound = 0; break; }
			}
			if (notfound)
			{
				term[maxterm] = j;
				cnt[maxterm] = 1;
				maxterm ++;
			}

		}

		double entropy = 0.0;
		for (k=0; k < maxterm; k++)
		{
			double p_k = ((double) cnt[k])/((double) nterms);
			entropy += p_k * log(p_k);
		}
		
		double y = -entropy;
		printf("%5d	%8.6g	%8.6g\n", i,x,y);

	}
}


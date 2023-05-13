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

	int total = 0;
	int totlen = 0;
	
	int p, q;
	for (q=2; q < 129; q++)
	{
		for (p=1; p<q; p++)
		{
			if (1 != gcf32(p,q)) continue;	
			f.SetRatio (p, q);

			int term[100];
			int cnt[100];
			int maxterm = 0;
			int nterms = f.GetNumTerms();
			totlen += nterms;
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
				double f_k = ((double) cnt[k])/((double) nterms);
				double p_k = term[k] + 1.0;
				p_k = -log(1.0 - 1.0 /(p_k * p_k));
				p_k /= M_LN2;
				entropy -= (f_k - p_k) * log(p_k);
			}
		
			double x = ((double) p)/ ((double) q);
			double y = entropy;
			// printf("%5d	%8.6g	%8.6g\n", p,x,y);
			printf("%8.6f	%5d	%5d	%8.6g\n", x, p,q,y);
			total ++;
		}

		double avglen = totlen / total;
		double ent = log(total) / M_LN2;
		ent /= avglen;
		printf("# total=%d avglen=%g ent=%g\n", total, avglen, ent);
	}
}


/*
 * irred-viz.c
 *
 * Find integer sequence for the golden polynomials.
 * And visualize it.  ... Nothing interesting happened here, dead end.
 * See irred-gold.c for basic utilities.
 *
 * February 2018, October 2020
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "irred-gold.c"

int main(int argc, char* argv[])
{
	int nmax = (1<<21);

	malloc_gold(nmax+1);

	int cnt = 0;
	int allcnt = 0;
	long psisum = 0;
	int plen = 0;
	int necklace = 1;

	int tord = 1<<plen;

	for (int n=0; n<nmax+1; n ++)
	{
		double gold = find_gold(n);

		// printf("---------\ngold=%g\n", gold);
		if (plen != order(n))
		{
			printf("# Total number for valid orbits of len=%d is %d\n", plen, cnt);
			plen = order(n);
			tord = 1<<plen;
			tord /= 2;
			cnt = 0;

			// advance and count necklace number
			necklace = 0;
			for (int k=n+1; k<nmax; k++)
			{
				double g = find_gold(k);
				if (0.5 < g) necklace++;
				if (order(k) != plen) break;
			}
			if (0 == necklace) break;
			printf("# Next necklace is %d\n", necklace);
		}
		if (0.5 < gold)
		{
			cnt ++;
			allcnt ++;
			psisum += n;

			// Always the case that tord/2 <= n < tord
			double frac = ((double) n) / ((double) tord);
			frac = 2.0*(frac - 0.5); // rescale to run 0 to 1.

			// OK, this should work.
			double seqfrac = ((double) cnt) / ((double) necklace);

			printf("%d	%d %d %d %g %g	%ld	%20.18g\n",
				allcnt, cnt, n, tord, seqfrac, frac, psisum, gold);
		}
	}
}

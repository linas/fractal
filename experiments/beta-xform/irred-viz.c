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

			// Always the case that tord/2 <= n < tord
			// double frac = ((double) n) / ((double) tord);
			// double frac = ((double) n) / ((double) (tord-1));
			double frac = ((double) n-1) / ((double) (tord-2));
			double ordfrac = 2.0*(frac - 0.5); // rescale to run 0 to 1.

			// First one is always zero, last is always one.
			// So e.g. necklace is 9, this will give eight steps
			double seqfrac = ((double) cnt-1) / ((double) necklace-1);
			// double seqfrac = ((double) cnt) / ((double) necklace);

			double ska = ((double) allcnt) /  ((double) tord);

			printf("%d	%d %d %d %g %g	%g	%20.18g\n",
				allcnt, cnt, n, tord, seqfrac, ordfrac, ska, gold);
		}
	}
}

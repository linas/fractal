/*
 * irred.c
 *
 * Find and verify irreducible golden polynomials.
 * February 2018
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "irred-gold.c"


/**
 * Find the zeros, then arrange them in sequential order.
 */
int bubble_sort(int nmax)
{
	// Find all of them.
	bool bad[nmax];
	for (int n=0; n<nmax; n ++)
	{
		double gold = find_gold(n);
		bad[n] = false;
		if (gold < 0.5) bad[n] = true;
	}

	// Discard the bad ones.
	int cnt = 0;

	// XXX Also discard the early ones.
	// for (int i=0; i< nmax; i++)
	for (int i=nmax/2; i< nmax; i++)
	{
		zero[cnt] = zero[i];
		if (!bad[i]) cnt++;
	}

	// Bubble sort them.
	for (int i=0; i< cnt; i++)
	{
		for (int j=i+1; j<cnt; j++)
		{
			if (zero[j] < zero[i])
			{
				double tmp = zero[j];
				zero[j] = zero[i];
				zero[i] = tmp;
			}
		}
	}

	return cnt;
}

#ifdef BIN_COUNT
int main(int argc, char* argv[])
{
	// int nmax = (1<<28) + 3;
	int nmax = (1<<18) + 3;
	// int nmax = (1<<20) + 3;

	setup_gold(nmax);

#define NPTS 1303
	int npts = NPTS;
	double bincnt[npts+1];
	for (int i=0; i<npts; i++)
	{
		bincnt[i] = 0.0;
	}

	int cnt = 0;
	for (int n=0; n<nmax; n ++)
	{
		double gold = find_gold(n);

		// Bin-count.
		if (gold < 0.5) continue;
		cnt ++;
		int nbin = (gold - 1.0) * npts;
		bincnt[nbin] += 1.0;
	}

	double norm = ((double) npts) / ((double) cnt);
	for (int i=0; i<npts; i++)
	{
		double x = ((double) i + 0.5) / ((double) npts);
		x += 1.0;
		double y = norm * bincnt[i];
		printf("%d	%g	%g\n", i, x, y);
	}
}

#endif // BIN_COUNT

int main(int argc, char* argv[])
{
	// int nmax = (1<<28) + 3;
	// int nmax = (1<<18) + 3;
	// int nmax = (1<<20) + 3;
	// int nmax = (1<<8) + 3;
	// int nmax = (1<<12) + 1;
	int nmax = (1<<20) + 1;

	malloc_gold(nmax);

// find_gold(16);
// printf("---------\ngold=%20.16g\n", zero[16]);
// exit(0);

#define PRINT_STUFF
#ifdef PRINT_STUFF
	int cnt = 0;
	int lyn = 1;
	int plen = 0;
	for (int n=0; n<nmax; n ++)
	{
		double gold = find_gold(n);

		// printf("---------\ngold=%g\n", gold);
		if (plen != order(n))
		{
			double rat = ((double) cnt) / ((double) lyn);
			printf("# total for len=%d is %d	ratio=%g\n", plen, cnt, rat);
			plen = order(n);
			lyn = cnt;
			cnt = 0;
		}
		if (0.5 < gold) { cnt++; }

		if (n < 128)
		// if (n < 540)
		// if (0)
		{
#if 0
			printf("%d l=%d %d %20.18g\n",
				n, len(n), cnt, gold);
#endif
			if (0.5 < gold)
			{
				printf("%d	%d	%20.18g #", n, order(n), gold);
				print_orbit(order(n), gold);
			}
		}
	}
#endif


#ifdef SORTED
	printf("#\n# Columns: idx, gold, ratio, ratio-bump\n#\n");
	int cnt = bubble_sort(nmax);
	for (int n=0; n<cnt-1; n ++)
	{
		double gold = zero[n];
		double ratio = zero[n+1] / gold;
		double bump = 1.0 / (ratio - 1.0);
		printf("%d	%g %g	%g\n", n, gold, ratio, bump);
	}
#endif
}

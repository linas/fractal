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

/* Return polynomial for bitstr */
double beta(int n, double x)
{
	double acc = 0.0;
	double xn = 1.0;
	int bitstr = 2*n+1;
	while (bitstr)
	{
		if (bitstr%2 == 1) acc += xn;
		xn *= x;
		bitstr >>= 1;
	}
	return xn - acc;
}

double find_zero(int n, double lo, double hi)
{
	double mid = 0.5 * (lo+hi);
	if (1.0e-15 > hi-lo) return mid;
	double fmid = beta(n, mid);
	if (0.0 < fmid) return find_zero(n, lo, mid);
	return find_zero(n, mid, hi);
}

/* Return length of bitstr, length in bits */
int len(int n)
{
	int len=0;
	int bitstr = 2*n+1;
	while (bitstr) { len++; bitstr >>= 1; }
	return len;
}

int main(int argc, char* argv[])
{
	int nmax = (1<<28) + 3;

#define EPS 2.0e-15

	int cnt = 0;
	int plen = 0;
	double* zero = (double*) malloc(nmax*sizeof(double));
	for (int n=0; n<nmax; n ++)
	{
		double gold = find_zero(n, 1.0, 2.0);
		zero[n] = gold;

#if 0
		bool aok = true;
		for (int j=0; j< n; j++)
		{
			double z = beta(n, zero[j]);
			if (fabs(z) < 4.0e-14) { aok = false; break; }
		}
#endif

		// printf("---------\ngold=%g\n", gold);
		bool ork = true;
		int nhl = n;
		int nh = nhl >> 1;
		while (nh)
		{
			// printf("duuude n=%d nhl=%d nh=%d znh=%g go=%g comp=%d bad=%d\n",
			// n, nhl, nh, zero[nh], gold, 0 == nhl%2, zero[nh] <= gold);
			if (0 == nhl%2 && zero[nh] < gold+EPS) {ork = false; break; }
			nhl = nh;
			nh >>= 1;
		}

		if (plen != len(n))
		{
			printf("total for len=%d is %d\n", plen, cnt);
			plen = len(n);
			cnt = 0;
		}
		if (ork) { cnt++; }

		if (n < 128)
		{
			printf("%d ok=%d l=%d %d %20.18g\n",
				n, ork, len(n), cnt, gold);
		}
	}
}

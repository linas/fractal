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

/* Return true, if its iterable */
bool iterable(int n)
{
	if (n <= 4) return true;
	return ((n-1)%4 != 0) && (iterable(n/2) || iterable(n/4)) ;
}

int main(int argc, char* argv[])
{
	int nmax = 64;

	int cnt = 0;
	int plen = 0;
	double zero[nmax];
	for (int n=0; n<nmax; n ++)
	{
		double gold = find_zero(n, 1.0, 2.0);
		zero[n] = gold;

		bool ok = true;
		for (int j=0; j< n; j++)
		{
			double z = beta(n, zero[j]);
			if (fabs(z) < 4.0e-14) { ok = false; break; }
		}
		if (plen != len(n)) {plen = len(n); cnt = 0; printf("\n");}
		if (ok && iterable(n)) cnt++;
		printf("%d ok=%d it=%d l=%d %d %20.18g\n",
			n, ok, iterable(n), len(n), cnt, gold);
	}
}

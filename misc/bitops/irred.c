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
	for (int n=0; n<20; n ++)
	{
		double gold = find_zero(n, 1.0, 2.0);
		printf("%d	%d	%20.18g\n", n, len(n), gold);
	}
}

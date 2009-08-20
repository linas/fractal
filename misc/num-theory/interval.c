
/*
 * interval.c
 *
 * replot ont an interval
 *
 * Linas Vepstas August 2009
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "moebius.h"

/* Reverse the bit order of p, assuming its of length n. */
int reverse(int p, int n)
{
	int s;
	int r = 0;
	for (s=0; s<n; s++)
	{
		r <<= 1;
		if (p & 0x1) r |= 0x1;
		p >>= 1;
	}
	return r;
}


int main(int argc, char * argv[])
{
	int k;

	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s len\n", argv[0]);
		exit(1);
	}

	int n = atoi(argv[1]);
	int m = 1<<n;
	double acc = 0.0;
	for (k=1; k<m; k++)
	{
		int r = reverse(k,n);
		double x = ((double) k) / ((double) m);

		// The most intersting ones, from the integral point of view, 
		// are mobius_mu and liouville lambda, since the integrate in
		// interesting ways.
		// double y = divisor(k);
		double y = moebius_mu(k);
		// double y = liouville_omega(k);
		// double y = liouville_lambda(k);
		// double y = mangoldt_lambda(k);

		acc += y;
		printf("%d	%d	%f	%f	%f\n", k, r, x, y, acc);
	} 

	return 0;
}

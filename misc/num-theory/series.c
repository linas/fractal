
/*
 * series.c
 *
 * graph of maclaurin series of totient
 *
 * Linas Vepstas December 2004
 */

#include <stdio.h>
#include "gcf.h"
#include "totient.h"

long double totient_series (long double x)
{
	long double acc = 0.0;

	long double xp = 1.0;
	int n=1;
	while (1)
	{
		acc += xp * totient_phi (n);

		if (xp < 1.0e-19) break;
		xp *= x;
		n++;
	}

	return acc;
}

int main ()
{
	int i;

	int nmax = 212;

	for (i=1; i<nmax; i++)
	{
		long double x = ((double) i)/((double) nmax);

		long double y = totient_series (x);
		y *= (1.0L-x)*(1.0L-x);

		printf ("%d	%Lg	%26.18Lg\n", i, x, y);
	}
}

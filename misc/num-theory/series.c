
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
	int n;
	for (n=1; n<1000000; n++)
	{
		acc += xp * totient_phi (n);

		if (xp < 1.0e-16) break;
		xp *= x;
	}

	return acc;
}

int main ()
{
	int i;

	int nmax = 212;

	for (i=1; i<nmax-10; i++)
	{
		double x = ((double) i)/((double) nmax);

		double y = totient_phi (nume);

		printf ("%d	%g	%g\n", i, x, y);
	}
}

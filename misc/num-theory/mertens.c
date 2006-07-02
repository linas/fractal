
/*
 * mertens.c
 *
 * Graph of merten function.
 * Ordinary x-y plot; output is ascii list of x-y values
 *
 * Linas Vepstas July 2006
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "moebius.h"

long double mertnes (int n)
{
	long double acc = 0.0;

	int i=1;
	for (i=1; i<=n; i++)
	{
		acc += moebius_mu (i);
	}

	return acc;
}

int main ()
{
	int i;

	int nmax = 41;

	for (i=1; i<nmax; i++)
	{
		long double y = mertens (x);

		printf ("%d	%26.18Lg\n", i, y);
		fflush (stdout);
	}
}

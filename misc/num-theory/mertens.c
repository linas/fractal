
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

int mertens (int n)
{
	int acc = 0.0;

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

	int nmax = 4100;

	int mert = 0;
	for (i=1; i<nmax; i++)
	{
		// mert = mertens (i);
		mert += moebius_mu(i);

		printf ("%d	%d\n", i, mert);
		fflush (stdout);
	}
}

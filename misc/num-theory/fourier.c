
/*
 * fourier.c
 *
 * Fourier transforms over varius number theoretic funcions
 *
 * Linas Vepstas November 2008
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

	int nmax = 10000000;
	nmax = 10000;
	
	int mert = 0;
	int scale = 1;
	for (i=1; i<nmax; i++)
	{
		// mert = mertens (i);
		mert += moebius_mu(i);

	}
}

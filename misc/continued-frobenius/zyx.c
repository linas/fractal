/*
 * zyx.c
 *
 * Graph the crazy-binomial function
 *
 * Linas Vepstas December 2004
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// #include <gsl/gsl_sf_zeta.h>
#include "zetafn.h"


/* Basic sum for integers in the middle of the sum */
long double symy (long double x, int m, long double ess)
{
	int i;

	long double en = -x;
	long double pnk = 1.0L;
	long double acc = 0.0L;

	for (i=0; i<50; i++)
	{
		long double eye = i;
		// long double term = gsl_sf_zeta (eye+ess);
		// term -= harmonic (m, eye+ess);
		long double term = zeta_minus_harmonic (m, eye+ess);
		term *= pnk;
		term *= fbinomial (ess+eye-1.0L, i);
		acc += term;
		pnk *= en;
	}
	
	return acc;
}


int main (int argc, char *argv[])
{
	int i;

	if (3>argc) {
		fprintf (stderr, "Usage: %s <ess> <m>\n", argv[0]);
		exit (1);
	}
	double ess = atof (argv[1]);
	int m = atoi (argv[2]);

	int nmax = 25;

	for (i=1; i<nmax; i++)
	{
		double x = ((double) i)/((double) nmax);
		x *= 7.0;

		double acc;
		acc = symy (x, m, ess);
		printf ("%d	%g	%g\n", i, x, acc);
	}

	return 0;
}

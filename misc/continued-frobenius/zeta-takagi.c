
/* zeta-takagi.c
 *
 * Goal: numerically validate a sum over zeta's
 * explore the weirdness of the takagi-related 
 * sums for non-integer values of s ... 
 * 
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_sf_zeta.h>
#include "zetafn.h"

/* Return the n'th harmnic number */
long double harmonic (int n, long double ess)
{
	long double acc = 0.0L;
	int k;

	for (k = 1; k<= n; k++)
	{
		long double kay = k;
		acc += powl (kay, -ess);
	}
	return acc;
}

/* Return Hurwitz zeta equiv of harmonic */
long double harmonic_hurwitz (int n, long double x, long double ess)
{
	long double acc = 0.0L;
	int k;

	for (k = 1; k<= n; k++)
	{
		long double kay = k;
		kay += x;
		acc += powl (kay, -ess);
	}
	return acc;
}

/* Return zeta function minus n'th harmonic */
long double zeta_minus_harmonic (int n, long double ess)
{
	long double acc = 0.0L;
	int k;

	for (k = n+1; k<10000; k++)
	{
		long double kay = k;
		long double term = powl (kay, -ess);
		acc += term;
		if (term < 1.0e-16*acc) break;
	}
	return acc;
}

long double sanity (long double ess)
{
	double zeta = gsl_sf_zeta (ess);
	double expect = harmonic (60, ess);

	return zeta-expect;
}

/* Most basic sum, with unity in the sum */
long double summy (long double ess)
{
	int i;

	long double twok = 1.0L;
	long double acc = 0.0L;
	long double ts = powl (0.5L, ess);

	for (i=0; i<50; i++)
	{
		long double eye = i;
		// long double term = zetam1(i+s) +1.0;
		long double term = gsl_sf_zeta (eye+ess);
		term *= twok;
		term *= fbinomial (ess+eye-1.0L, i);
		acc += term;
		twok *= -0.5L;
	}
	acc *= ts;
	acc += 1.0;
	// double expect = zetam1(s) +1.0;
	double expect = gsl_sf_zeta (ess);
	expect *= 1.0-ts;
	acc -=  expect;
	
	return acc;
}

/* Basic sum for integers in the middle of the sum */
long double symy (int m, int n, long double ess)
{
	int i;

	long double en = -n;
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
	// double expect = gsl_sf_zeta (ess);
	// expect -= harmonic (m+n, ess);
	double expect = zeta_minus_harmonic (m+n, ess);
	acc -=  expect;
	
	return acc;
}

/* Basic half-integer sum */
long double sym_half (int m, int n, long double ess)
{
	int i;

	long double en = n;
	long double pnk = 1.0L;
	long double acc = 0.0L;
	long double ts = powl (0.5L, ess);

	for (i=0; i<60; i++)
	{
		long double eye = i;
		// long double term = gsl_sf_zeta (eye+ess);
		// term -= harmonic (m, eye+ess);
		long double term = zeta_minus_harmonic (m, eye+ess);
		term *= pnk;
		term *= fbinomial (ess+eye-1.0L, i);
		acc += term;
		pnk *= -en*(1.0 + 0.5*en);
	}
	acc *= ts;

	double expect = gsl_sf_zeta (ess);
	expect *= 1.0 - ts;
	expect -= 1.0;
	expect -= ts * harmonic_hurwitz (m+n, 0.5L, ess);
	acc -=  expect;
	
	return acc;
}

/* half-integer sum */
long double sym_half_gen (int m, int n, long double ess)
{
	int i;

	long double en = n;
	long double pnk = 1.0L;
	long double acc = 0.0L;
	long double ts = powl (0.5L, ess);

	for (i=0; i<70; i++)
	{
		long double eye = i;
		// long double term = gsl_sf_zeta (eye+ess);
		// term -= harmonic (m, eye+ess);
		long double term = zeta_minus_harmonic (m, eye+ess);
		term *= pnk;
		term *= fbinomial (ess+eye-1.0L, i);
		acc += term;
printf ("duude acc=%Lg  term=%Lg  nk*bin=%Lg\n", acc, term, pnk*fbinomial (ess+eye-1.0L, i));
		pnk *= -(en+0.5);
	}
	acc *= ts;

	double expect = gsl_sf_zeta (ess);
	expect *= 1.0-ts;
	expect -= harmonic (2*n+2*m+1, ess);
	expect += ts*harmonic (n+m, ess);
	acc -=  expect;
	
	return acc;
}

int
main (int argc, char * argv[])
{
	int i;

	if (3>argc) {
		fprintf (stderr, "Usage: %s <m> <n>\n", argv[0]);
		exit (1);
	}
	int m = atoi (argv[1]);
	int n = atoi (argv[2]);

	int nmax = 25;

	for (i=0; i<nmax; i++)
	{
		double x = ((double) i)/((double) nmax);
		x *= 6.0;
		x += 1.3;

		double acc;
		// acc = summy (x);
		// acc = symy (6, 7, x);
		// acc = sym_half (m, n, x);
		acc = sym_half_gen (m, n, x);
		printf ("%d	%g	%g\n", i, x, acc);
	}

	return 0;
}

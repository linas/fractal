
/* zeta-takagi.c
 *
 * Goal: numerically validate a sum over zeta's
 * explore the weirdness of the takagi-related 
 * sums for non-integer values of s ... 
 * 
 */

#include <math.h>
#include <stdio.h>

#include <gsl/gsl_sf_zeta.h>
#include "zetafn.h"

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

long double symy (int m, int n, long double ess)
{
	int i;

	long double en = -n;
	long double pnk = 1.0L;
	long double acc = 0.0L;

	for (i=0; i<30; i++)
	{
		long double eye = i;
		long double term = gsl_sf_zeta (eye+ess);
		term -= harmonic (m, eye+ess);
		term = zeta_minus_harmonic (m, eye+ess);
		term *= pnk;
		term *= fbinomial (ess+eye-1.0L, i);
		acc += term;
		pnk *= en;
	}
	double expect = gsl_sf_zeta (ess);
	expect -= harmonic (m+n, ess);
	expect = zeta_minus_harmonic (m+n, ess);
	acc -=  expect;
	
	return acc;
}

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

int
main (int argc, char * argv[])
{
	int i;

	int nmax = 25;

	for (i=0; i<nmax; i++)
	{
		double x = ((double) i)/((double) nmax);
		x *= 6.0;
		x += 1.3;

		double acc = summy (x);
		acc = symy (6, 7, x);
		printf ("%d	%g	%g\n", i, x, acc);
	}

	return 0;
}


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
long double sym_integer (int m, int n, long double ess)
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


/* Generic half-integer sum i.e. for n+1/2  */
long double sym_half (int m, int n, long double ess)
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
		pnk *= -(en+0.5);
	}
	acc *= ts;

	double expect = gsl_sf_zeta (ess);
	expect *= 1.0-ts;
	expect -= harmonic (2*n+2*m+1, ess);
	expect += ts*harmonic (n+m, ess);
	acc -=  expect;
	
/*
	This is an equivalent to above.
	double expect = gsl_sf_zeta (ess);
	expect *= 1.0 - ts;
	expect -= 1.0;
	expect -= ts * harmonic_hurwitz (m+n, 0.5L, ess);
	acc -=  expect;
*/
	return acc;
}

/* The 3/4ths reflective sum */
long double reflect_fourths (long double ess)
{
	int i;

	long double fk = 1.0L;
	long double tk = 1.0L;
	long double acc = 0.0L;

	for (i=0; i<70; i++)
	{
		long double eye = i;
		long double term = gsl_sf_zeta (eye+ess);
		term *= fk +tk;
		term *= fbinomial (ess+eye-1.0L, i);
		acc += term;
		fk *= -0.25;
		tk *= -0.75;
	}
	acc *= powl (4.0, -ess);
	acc += powl (3.0, -ess);
	acc += 1.0;

	double expect = gsl_sf_zeta (ess);
	expect *= 1.0 - pow (2.0, -ess);
	acc -=  expect;
	
	return acc;
}

/* The general 3/4ths sum */
long double sym_fourths (int n, long double ess)
{
	int i;

	long double en = 0.5*n;
	long double fk = 1.0L;
	long double tk = 1.0L;
	long double acc = 0.0L;
	long double ts = powl (2.0, -ess);

	for (i=0; i<70; i++)
	{
		long double eye = i;
		long double term = gsl_sf_zeta (eye+ess);
		term *= fk+tk;
		term *= fbinomial (ess+eye-1.0L, i);
		acc += term;
		fk *= -(en+0.25);
		tk *= -(en+0.75);
	}
	acc *= powl (4.0, -ess);
	long double alt = harmonic (2*n+3, ess);
	alt -= ts*harmonic (n+1, ess);

	acc += alt;

	double expect = gsl_sf_zeta (ess);
	expect *= 1.0-ts;
	acc -=  expect;
	
	return acc;
}


double check (int n, double ess)
{
	double ts = pow (2.0, -ess);
	double expect = 1.0 + ts * harmonic_hurwitz (n, 0.5L, ess);

	expect -= harmonic (2*n+1, ess);
	expect += ts*harmonic (n, ess);

	return expect;
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

		double acc, a2=0.0;
		// acc = summy (x);
		// acc = sym_integer (6, 7, x);
		// acc = sym_half (m, n, x);
		acc = sym_fourths (n, x);
		acc = reflect_fourths (x);
		printf ("%d	%g	%g\n", i, x, acc);
	}

	return 0;
}

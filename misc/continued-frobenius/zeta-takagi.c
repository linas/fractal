
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

long double summy (long double ess)
{
	int i;

	long double twok = 1.0L;
	long double acc = 0.0L;
	int s = ess;
	long double ts = powl (0.5L, ess);

	for (i=0; i<50; i++)
	{
		// long double term = zetam1(i+s) +1.0;
		long double term = gsl_sf_zeta (i+ess);
		term *= twok;
		term *= binomial (s+i-1, s-1);
		acc += term;
		twok *= -0.5;
	}
	acc *= ts;
	acc += 1.0;
	acc /= 1.0-ts;
	// acc -= zetam1(s) +1.0;
	acc -= gsl_sf_zeta (ess);
	
	return acc;
}

int
main ()
{
	int i;

	int nmax = 353;

	for (i=0; i<nmax; i++)
	{
		double x = ((double) i)/((double) nmax);
		x *= 6.0;
		x += 1.3;

		double acc = summy (x);
		printf ("%d	%g	%g\n", i, x, acc);
	}

	return 0;
}

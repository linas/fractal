/*
 * exact-series.c
 *
 * Solution via series
 *
 * Linas Vepstas December 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "binomial.h"

/* Integral of (x^(s-1))/(x+alpha) dx
 */
long double complex eff (long double complex s, long double alpha, long double x)
{
	int k;
	long double complex term, sum;
	long double opxa = 1.0 + x/alpha;
	long double opxak = opxa;

printf("oxpa = %Lg\n", opxa);

	sum = 0.0L;

	for (k=1; k<155456123; k++)
	{
		term = cbinomial(s,k);
		term *= opxak;
		term /= (long double) k;
		if (k%2)
		{
			sum -= term;
		}
		else
		{
			sum += term;
		}

		double tm = cabs(term);
		if (tm < 1.0e-20) break;

		opxak *= opxa;
	}

	long double xa = x+alpha;
	if (xa <= 0.0L)
	{
		fprintf(stderr, "Error: unexpected sign for x+a=%Lg\n", xa);
		exit(1);
	}
	term = logl(-xa);
	sum += term;

	term = s * logl(-alpha);
	term = cexpl(term);

	sum *= term;

	return sum;
}

/* Integral of s_12 */

long double complex gral_s12(long double complex s, unsigned int a1max, unsigned int a2max)
{
	unsigned int na1, na2;
	long double complex term;
	long double complex sum = 0.0L;
	
	for (na1=1; na1<a1max; na1++)
	{
		for (na2=1; na2<a2max; na2++)
		{
			long double a1 = na1;
			long double a2 = na2;
			long double b = a1 - a2;
			long double c = - (a1 * a2 + 1.0L) * b;
			long double d = 1.0L + a2 * b;
			long double a = 1.0L - a1 * b;
			long double xlo = a2 / (1.0L + a1 * a2);
			long double xhi = (1.0L + a2) / (1.0L + a1 + a1 * a2);
			long double greb = d / c;

			term = eff(s, greb, xlo);
			sum += term;
			term = eff(s, greb, xhi);
			sum += term;
		}
	}
	return sum;
}

main ()
{
	long double complex ess = 0.5 + I*12.0;

	long double complex ans = gral_s12(ess, 100, 100);
printf("ans= %g %g\n\n", creal(ans), cimag(ans));

}

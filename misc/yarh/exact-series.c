/*
 * exact-series.c
 *
 * Solution via series
 *
 * Linas Vepstas December 2010
 */

#include <stdio.h>
#include <math.h>

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

	sum = 0.0;

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

long double complex sum_sum(long complex s)
{

}

main ()
{

	int na1 = 2;
	int na2 = 1;

			long double a1 = na1;
			long double a2 = na2;
			long double b = a1 - a2;
			long double c = - (a1 * a2 + 1.0L) * b;
			long double d = 1.0L + a2 * b;
			long double a = 1.0L - a1 * b;
			long double xlo = a2 / (1.0L + a1 * a2);
			long double xhi = (1.0L + a2) / (1.0L + a1 + a1 * a2);
			long double greb = d / c;

	printf("duude a1=%d a2=%d b=%Lg c=%Lg d=%Lg\n", na1, na2, b,c,d);

	long double complex ess = 0.5 + I*12.0;

printf("duude xlo=%Lg xhi=%Lg alpha=%Lg\n", xlo, xhi, greb);
	
	long double complex ans = eff(ess, greb, xlo);
printf("ans= %g %g\n\n", creal(ans), cimag(ans));
	ans = eff(ess, greb, xhi);
printf("ans= %g %g\n", creal(ans), cimag(ans));

}

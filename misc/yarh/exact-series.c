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
	long double opxa = 1.0+x/alpha;
	long double opxak = opxa;

	sum = 0.0;

	for (k=1; k<100; k++)
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

printf("duuude k=%d re=%f im=%f\n", k, creal(sum), cimag(sum));

		opxak *= opxa;
	}

	return sum;
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
			long double xlo = a2 / (1.0L + a1 * a2);
			long double xhi = (1.0L + a2) / (1.0L + a1 + a1 * a2);
			long double greb = d / c;

	printf("duude a1=%d a2=%d b=%Lg c=%Lg d=%Lg\n", na1, na2, b,c,d);

	long complex ess = 0.5 + I*12.0;

printf("duude xlo=%Lg xhi=%Lg alpha=%Lg\n", xlo, xhi, greb);
	
	eff(ess, greb, xlo);
	eff(ess, greb, xhi);

}

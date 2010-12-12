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

printf("duuude k=%d re=%llf im=%llf\n", k, creal(sum), cimag(sum));

		opxak *= opxa;
	}

	return sum;
}

main ()
{

	int a1 = 1;
	int a2 = 2

	double a = 1+a1*(a2-a1);
	double b = a1-a2;
	double c = (a1*a2+1)*(a2-a1);
	double d = 1+a2*(a1-a2);

	double alp = d/c;

	long complex ess = 0.5 + I*12.0;

	long double x = (1+a2)/(1+a1+a1*a2);

pritf("duude x=%g alpha=%g\n", x, alp);
	
	eff(ess, alp, x);

}

/*
 * newton-series.c
 *
 * Newton-series summation of the continued fraction
 * riemann integral. For S_12 only.
 *
 * Linas Vepstas December 2010
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

complex alpha(complex ess, int a1max, int a2max) 
{
	unsigned int na1, na2;
	complex sum = 0.0;

	for (na1=1; na1<a1max; na1++)
	{
		for (na2=1; na2<a2max; na2++)
		{
			if (na1 == na2) continue;  // c==0, d==0

			double a1 = na1;
			double a2 = na2;
			double c = (a1 * a2 + 1.0) * (a2 - a1);
			double d = 1.0 + a2 * (a1 - a2);
			double xlo = a2 / (1.0 + a1 * a2);
			double xhi = (1.0 + a2) / (1.0 + a1 + a1 * a2);
			complex beta = - ((double) d) / ((double) c);

			complex term = log((xhi+beta) /(xlo+beta));
			term *= cpow(beta, ess-1.0);
			term /= c*c ;
			sum += term;
		}
	}
	return sum;
}

int main (int argc, char * argv[])
{

	complex ess = 0.5 + I*7.0;
	complex yo = alpha(ess, 121, 121);

	double re = creal (yo);
	double im = cimag(yo);

	printf("duude %g %g\n", re, im);
	return 0;
}

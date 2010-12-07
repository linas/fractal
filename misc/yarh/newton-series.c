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
	unsigned int a1, a2;
	complex sum = 0.0;

	for (a1=1; a1<a1max; a1++)
	{
		for (a2=1; a2<a2max; a2++)
		{
			int c = (a1 * a2 + 1) * (a2 - a1);
			int d = 1 + a2 * (a1 - a2);
			complex beta = - ((double) d) / ((double) c);
			complex term = cpow(beta, ess-1.0);
			sum += term;
		}
	}
	return sum;
}

int main (int argc, char * argv[])
{

	complex ess = 0.5 + 7.0;
	complex yo = alpha(ess, 1123, 1123);

	double re = creal (yo);
	double im = cimag(yo);

	printf("duude %g %g\n", re, im);
	return 0;
}


/*
 * series.c
 *
 * graph of maclaurin series of totient
 *
 * Linas Vepstas December 2004
 */

#include <math.h>
#include <stdio.h>

#include "divisor.h"
#include "gcf.h"
#include "totient.h"

long double totient_series (long double x)
{
	long double acc = 0.0;

	long double xp = 1.0;
	int n=1;
	while (1)
	{
		long double term = xp * totient_phi (n);
		acc += term;

		if (term < 1.0e-20*acc) break;
		xp *= x;
		n++;
	}

	return acc;
}

// return the limit as the totient sum goes to x-> 1
void limit (void)
{
	long double p = 0.5L;
	long double prev = 0.0;
	int i=1;
	while (1)
	{
		long double x = 1.0L - p;

		long double y = totient_series (x);
		y *= p*p;

		long double guess = y + (y-prev)/3.0L;
		printf ("%d	%Lg	%26.18Lg	%26.18Lg\n", i, x, y,  guess);

		p *= 0.5L;
		i++;
		prev = y;
	}
}

long double divisor_series (long double x)
{
	long double acc = 0.0;

	long double xp = 1.0;
	int n=1;
	while (1)
	{
		long double term = xp * divisor (n);
		acc += term;

		if (term < 1.0e-20*acc) break;
		xp *= x;
		n++;
	}

	return acc;
}

int main ()
{
	int i;

	int nmax = 432;

	for (i=1; i<nmax; i++)
	{
		long double x = ((double) i)/((double) nmax);

#define TOTIENT_SERIES
#ifdef TOTIENT_SERIES
		long double y = totient_series (x);
		y *= (1.0L-x)*(1.0L-x);

		long double z = 0.607927101 * sin (0.5*M_PI*x);
		printf ("%d	%Lg	%26.18Lg	%26.18Lg\n", i, x, y, z);
#endif


#if 0
		long double y = divisor_series (x);
		// y *= (1.0L-x)*(1.0L-x);
		printf ("%d	%Lg	%26.18Lg\n", i, x, y);
#endif

	}
}


/*
 * angular.c
 *
 * Graphs of maclaurin series of totient and other number-theoretic
 * arithmetic series.  These are plots taken along a circular slice;
 * output is ascii list of theta-f(theta) values
 *
 * Linas Vepstas December 2004
 * Linas Vepstas May 2005
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "divisor.h"
#include "gcf.h"
#include "moebius.h"
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

long double c_divisor_series (long double x)
{
	long double complex acc = 0.0;

	long double complex xi = x*I;
	long double complex xp = x*I;
	int n=1;
	while (1)
	{
		long double complex term = xp * divisor (n);
		acc += term;

		if (cabsl(term) < 1.0e-20*cabsl(acc)) break;
		xp *= xi;
		n++;
	}

	return cabsl (acc);
}

long double c_erdos_series (long double x)
{
	long double complex acc = 0.0;

	long double complex xi = x*I;
	long double complex xp = x*I;
	while (1)
	{
		long double complex term = xp / (1.0L-xp);
		acc += term;

		if (cabsl(term) < 1.0e-20*cabsl(acc)) break;
		xp *= xi;
	}

	return cabsl(acc);
}

long double erdos_series (long double x)
{
	long double acc = 0.0;

	long double xp = x;
	while (1)
	{
		// long double term = xp / (1.0L-xp);
		long double term = 1.0L/(xp * (xp-1.0L));
		acc += term;

		if (term < 1.0e-20*acc) break;
		xp *= x;
	}

	return acc;
}

long double z_erdos_series (long double x)
{
	long double acc = 0.0;

	long double tk = 0.5L;
	long double xp = x;
	while (1)
	{
		long double term = xp / (1-tk);
		acc += term;

		if (term < 1.0e-20*acc) break;
		xp *= x;
		tk *= 0.5L;
	}

	return acc;
}

long double moebius_series (long double x)
{
	long double acc = 0.0;

	long double xp = 1.0;
	int n=1;
	while (1)
	{
		long double term = xp * moebius_mu (n);
		acc += term;

		if (xp < 1.0e-18) break;
		xp *= x;
		n++;
	}

	return acc;
}

int main ()
{
	int i;

	int nmax = 41;

	long double tp = 0.5;
	for (i=1; i<nmax; i++)
	{
		long double x = ((double) i)/((double) nmax);

// #define TOTIENT_SERIES
#ifdef TOTIENT_SERIES
		long double y = totient_series (x);
		y *= (1.0L-x)*(1.0L-x);
		long double z = 0.607927101 * sin (0.5*M_PI*x);

		long double r = 2.0L*(y/z - 1.0L);
		printf ("%d	%Lg	%26.18Lg	%26.18Lg	%26.18Lg\n", i, x, y, z,r);
#endif


// #define DIVISOR_SERIES
#ifdef DIVISOR_SERIES
		x = 1.0L-tp;
		long double y = divisor_series (x);
		// y *= (1.0L-x)*(1.0L-x);
		y *= tp;

		printf ("%d	%Lg	%26.18Lg\n", i, x, y);
#endif

// #define C_DIVISOR_SERIES
#ifdef C_DIVISOR_SERIES
		long double y = c_divisor_series (x);
		long double z = c_erdos_series (x);

		printf ("%d	%Lg	%26.18Lg	%26.18Lg	%26.18Lg\n", i, x, y, z, y-z);
#endif

// #define ERDOS_SERIES
#ifdef ERDOS_SERIES
		long double y = z_erdos_series (x);
		y *= (1.0L-x);

		printf ("%d	%Lg	%26.18Lg\n", i, x, y);
#endif

// #define MOEBIUS_SERIES
#ifdef MOEBIUS_SERIES
		long double y = moebius_series (x);

		printf ("%d	%Lg	%26.18Lg\n", i, x, y);
		fflush (stdout);
#endif


		tp *= 0.5L;
	}
}

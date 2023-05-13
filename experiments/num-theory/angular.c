
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
#include <stdlib.h>

#include "gcf.h"
#include "moebius.h"
#include "modular.h"
#include "totient.h"


long double divisor_series (long double r, long double theta)
{
	long double re_acc = 0.0;
	long double im_acc = 0.0;

	long double rp = 1.0;
	long double thp = 0.0;
	int n=1;
	while (1)
	{
		// long double term = rp * divisor (n);
		long double term = rp * sigma (n,2);
		// long double term = rp * totient_phi (n);
		// long double term = xp * moebius_mu (n);
		long double re_term = term * cosl (thp);
		long double im_term = term * sinl (thp);
		re_acc += re_term;
		im_acc += im_term;

		if (rp < 1.0e-20) break;
		rp *= r;
		thp += theta;
		n++;
	}

	long double acc = sqrtl (re_acc*re_acc + im_acc*im_acc);
	// return re_acc;
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


int 
main (int argc, char * argv[])
{
	int i;

	if (argc < 3) {
		fprintf (stderr, "Usage: %s <1/(1-radius)> <npts>\n", argv[0]);
		exit (1);
	}
	
	long double radius = atof(argv[1]);
	radius = 1.0L - 1.0L/radius;
	int nmax = atoi (argv[2]);

	printf ("#\n# angular series\n#\n");
	printf ("# radius=%Lg  npts=%d\n#\n", radius, nmax);

	long double acc = 0.0;
	for (i=1; i<nmax; i++)
	{
		long double x = ((double) i)/((double) nmax);

		long double y = divisor_series (radius, 2.0*M_PI* x);
		acc += 1.0L/y;

		printf ("%d	%Lg	%26.18Lg	%Lg\n", i, x, y, acc);

	}
}

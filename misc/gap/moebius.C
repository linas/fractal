
/* 
 * moebius.C
 *
 * moebius transform of the Takagi curve
 *
 * Linas October 2004
 */

#include <complex.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

long double triangle (long double x)
{
	long double t = x - floorl(x);
	if (0.5L > t) return 2.0L*t;
	return 2.0L*(1.0L-t);
}


// the farey/isola map
long double pointy (long double x)
{
	long double t = x - floorl(x);
	if (0.5L < t) return (1.0L-t)/t;
	return t/(1.0L-t);
}

long double arith_triangle (long double x, int n)
{
	long double kay = n-1;
	long double tp = powl (2.0, kay);
	// long double fn = triangle (tp*x);
	long double fn = pointy (tp*x);
	return fn;
}

long double moebius (long double x, int n)
{
	int d;

	long double acc = 0.0;
	for (d=1; d<= n; d++)
	{
		if (n%d) continue;
		acc += arith_triangle (x, d);
	}
	return acc;
}


/* The moebius-transformed takagi curve */
long double moebius_takagi (long double w, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	for (k=1; k<50; k++)
	{
		long double term = tw* moebius (x, k);
		acc += term;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	acc *= 1.0-w;
	return acc;
}

/* The dirichlet-summed takagi curve */
long double moebius_dirichlet (long double w, long double x)
{
	int k;
	long double acc = 0.0L;
	for (k=1; k<50; k++)
	{
		long double kay = k;
		long double tw = powl (kay, -w);
		long double term = tw*moebius (x, k);
		acc += term;
		if (1.0e-16 > tw) break;
	}

	return acc;
}


main (int argc, char *argv[])
{
	int i;

	int nmax = 2513;

	if (argc <2)
	{
		printf ("Usage: %s <w-value>\n", argv[0]);
		exit (1);
	}
	double w = atof(argv[1]);

	for (i=0; i<nmax; i++)
	{
		double x = i/((double)nmax);
		double tw = moebius_takagi (w, x);
		// double tw = moebius_dirichlet (w, x);


		printf ("%d	%8.6g	%8.6g\n", i, x, tw);
		fflush (stdout);
	}
}

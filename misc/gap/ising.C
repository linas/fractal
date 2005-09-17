
/* 
 * ising.C
 *
 * Ising model inspired takagi curve stuff 
 *
 * Linas October 2004
 * Linas September 2005
 */

#include <complex.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Farey.h"
ContinuedFraction far;

long double triangle (long double x)
{
	long double t = x - floorl(x);
	if (0.5L > t) return 2.0L*t;
	return 2.0L*(1.0L-t);
}

long double balanced_triangle (long double x)
{
	return (triangle(x) - 0.5L);
}

/* The main, core basic takagi curve */
long double takagi (long double w, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	long double tp = 1.0L;
	for (k=0; k<50; k++)
	{
		long double term = tw* balanced_triangle (tp*x);
		acc += term;
		tp *= 2.0L;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

// ===============================================
// ===============================================

main (int argc, char *argv[])
{
	int i;

	// int nmax = 512;
	// int nmax = 432;
	// int nmax = 1717;
	int nmax = 2048;

	if (argc <2)
	{
		printf ("Usage: %s <w-value>\n", argv[0]);
		exit (1);
	}
	double w = atof(argv[1]);

	double acc = 0.0;
	for (i=0; i<nmax; i++)
	{
		double x = i/((double)nmax);

		/* jitter, as this can make a difference */
		// x += ((double) rand()) / (RAND_MAX*((double)nmax));
		
		// double ts = isola (w, x);
		double tw = takagi (w, x);
		double ts = tw;
		tw = exp (-tw);
		acc += tw;
		ts = acc;

		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, tw, ts, tw-ts);
		fflush (stdout);
	}
}

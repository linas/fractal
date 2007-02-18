/* 
 * takagi-n.C
 *
 * Draw non-dyadic variant of the Takagi curve
 *
 * Linas December 2005
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

long double balanced_triangle (long double x)
{
	return (triangle(x) - 0.5L);
}

double square (double x)
{
	x -= floor (x);
	if (1.0 > x*2.0) return -2.0;
	return 2.0;
}

long double parabola_down (long double x)
{
	long double t = x - floorl(x);
	return 4.0L*t*(1.0L-t);
}

long double parabola_up (long double x)
{
	long double t = x - floorl(x);
	if (0.5 < t) t = 1.0-t;
	return 4.0L*t*t;
}

long double sine (long double x)
{
	long double t = sinl (M_PI*x);
	return t*t;
}

/* Have triangle bumps at all integer frequencies */
long double n_takagi (long double w, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	for (k=0; k<350; k++)
	{
		int d;
		long double term = tw* triangle ((k+1)*x);
		// long double term = tw* sine ((k+1)*x);
		acc += term;
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
	int nmax = 721;
	// int nmax = 432;
	// int nmax = 1717;
	// int nmax = 2048;

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
		
		double tw = n_takagi (w, x);
		tw *= (1.0-w);
		tw = 1.0-tw;
		double ts = log(tw);
		ts = tw;

		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, tw, ts, tw-ts);
		fflush (stdout);
	}
}

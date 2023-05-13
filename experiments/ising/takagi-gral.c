/* 
 * takagi-gral.c
 *
 * draw the integrals of the Takagi curve
 *
 * Linas October 2004
 * Linas September 2005
 */

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

/* The main, core basic takagi curve */
long double takagi (long double w, long double x)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	long double tp = 1.0L;
	for (k=0; k<50; k++)
	{
		long double term = tw* triangle (tp*x);
		// long double term = tw* balanced_triangle (tp*x);
		// long double term = tw* parabola_down (tp*x);
		// long double term = tw* parabola_up (tp*x);
		acc += term;
		tp *= 2.0L;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

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

	double delta = 1.0/((double)nmax);
	double acc = 0.0;
	double acc2 = 0.0;
	for (i=0; i<nmax; i++)
	{
		double x = i*delta;

		/* jitter, as this can make a difference */
		// x += ((double) rand()) / (RAND_MAX*((double)nmax));
		
		double tw = takagi (w, x);
		double ts = exp (-tw);
		acc += ts*delta;

		/* and now, by half .. */ 
		if (i%2 == 0) {
			tw = takagi (w, 0.5*x);
			ts = exp (-tw);
			acc2 += ts*delta;
		}

		double fit = acc2;
		
		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, acc, fit, fit-acc);
		fflush (stdout);
	}
}

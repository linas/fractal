
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
		// long double term = tw* balanced_triangle (tp*x);
		long double term = tw* triangle (tp*x);
		acc += term;
		tp *= 2.0L;
		tw *= w;
		if (1.0e-16 > tw) break;
	}

	return acc;
}

/* x = p / 2**n */
int itriangle (int p, int n)
{
	unsigned int deno = 1<<n;

	if (p >= deno) p-= deno;
	if (p == 0) return 0;

	if (2*p < deno) return p;

	p = deno - p;
	return p;
}

/* x = p / 2**n */
long double itakagi (long double w, int p, int n)
{
	int k;
	long double acc = 0.0L;
	long double tw = 1.0L;
	long double tp = 1.0L;
	for (k=0; k<n; k++)
	{
		p = itriangle (p, n-k);
		if (0 == p) break;
		long double term = tw * ((long double) p) / ((long double) (1<<(n-k)));
		acc += term;
		tp *= 2.0L;
		tw *= w;
	}

	return acc;
}

// ===============================================
// ===============================================

main (int argc, char *argv[])
{
	int i;

	if (argc <2)
	{
		printf ("Usage: %s <len> <w-value>\n", argv[0]);
		exit (1);
	}
	int p = atoi (argv[1]);
	double w = atof(argv[2]);
	
	// int nmax = 512;
	// int nmax = 432;
	// int nmax = 1717;
	// int nmax = 2048;
	// int nmax = 8192;
	// int nmax = 32768;
	// int nmax = 32;

	int nmax = 1<<p;
	
	double scale = sqrt ((double) nmax);
	printf ("#\n# scale=%g\n#\n", scale);

	double fa = 0.0;
	double acc = 0.0;
	double xlast = 0.0;
	double step = 0.005;
	for (i=1; i<nmax; i++)
	{
		double x = i/((double)nmax);

		/* jitter, as this can make a difference */
		// x += ((double) rand()) / (RAND_MAX*((double)nmax));
		
		// double ts = isola (w, x);
		// double tw = takagi (w, x);
		double tw = itakagi (w, i, p);
		
		// ts = exp (-2.0*0.69314*tw);
		double ts = pow (4.0, -tw);
		// ts = pow (M_PI, -tw);
		// ts = pow (9.0/4.0, -tw);
		acc += ts;

		if (x> xlast + step) {
			fa = InvFarey(x);
			xlast += step;
		}

		printf ("%d	%8.6g	%8.6g	%8.6g	%8.6g\n", i, x, tw, acc, fa);
		fflush (stdout);
	}
}

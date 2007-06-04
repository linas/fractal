
/*
 * liouville.c
 *
 * Graphs of series involving Liouville function.
 * These are ordinary x-y plots; output is 
 * ascii list of x-y values
 *
 * Linas Vepstas July 2006
 */

#include <math.h>
#include <stdio.h>

#include "modular.h"

int main ()
{
	int i;

	/* polya conjecture broken at n=906150257 also n=906180359 */
	// int nmax = 10000;
	unsigned long nmax = (1<<31) + (1<<30) + (1<<29);
	double logscale = 2.0;
	int scale = logscale;

	/* npoints = number of points in the graph */
	int npoints = 900;
	double scale_delta = exp (log((double)nmax) / npoints);
	printf ("#\n# step scale factor = %g\n#\n", scale_delta);

	int sum = 0;
	long double lead=0.0L;
	long double inf = 1.0e100;
	long double sup = -1.0e100;
	for (i=1; i<=nmax; i++)
	// for (i=nmax; i>0; i--)
	{
		long double x = ((double) i)/((double) nmax);

		int d = liouville_lambda (i);
		// int d = liouville_omega (i);
		sum += d;
		// lead += d / ((long double) i);
		lead = sum;

		if (inf>lead) inf=lead;
		if (sup<lead) sup=lead;
		if (0 == i%scale)
		{
			// lead = sum;
			// lead -= 0.5*i*logl(i);
			// printf ("%d	%d	%26.18Lg\n", i, d, lead);

			/* print locations where conjecture is false */
			printf ("%d	%d	%26.18Lg\n", i, d, inf);
			if (0 < sup)
			{
				printf ("%d	%d	-1\n", i, d);
				printf ("%d	%d	1\n", i, d);
			}
			printf ("%d	%d	%26.18Lg\n", i, d, sup);
			fflush (stdout);
			inf = 1.0e100;
			sup = -1.0e100;

			logscale *= scale_delta;
			scale = logscale;
		}
	}
}

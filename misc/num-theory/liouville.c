
/*
 * liouville.c
 *
 * Graphs of series involving lioville function.
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

	int nmax = 1000;
	int scale = 1;

	int sum = 0;
	long double lead=0.0L;
	for (i=1; i<=nmax; i++)
	// for (i=nmax; i>0; i--)
	{
		long double x = ((double) i)/((double) nmax);

		int d = liouville_lambda (i);
		// int d = liouville_omega (i);
		sum += d;
		lead += d / ((long double) i);
		if (0 == i%scale)
		{
			// lead = sum;
			// lead -= 0.5*i*logl(i);
			printf ("%d	%d	%26.18Lg\n", i, d, lead);
		}
	}
}

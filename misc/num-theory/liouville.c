
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

	int nmax = 2500;
	int scale = 1;

	int sum = 0;
	for (i=1; i<=nmax; i++)
	// for (i=nmax; i>0; i--)
	{
		long double x = ((double) i)/((double) nmax);

		int d = liouville_lambda (i);
		sum += d;
		if (0 == i%scale)
		{
			long double lead = sum;
			printf ("%d	%d	%26.18Lg\n", i, d, lead);
		}
	}
}

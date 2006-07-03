
/*
 * von_mangoldt.c
 *
 * Graph of Chebyshev (sum of von Mangoldt) function.
 * Ordinary x-y plot; output is ascii list of x-y values
 *
 * Linas Vepstas July 2006
 */

#include <complex.h>
#include <math.h>
#include <stdio.h>

#include "moebius.h"

long double chebyshev (int n)
{
	long double acc = 0.0;

	int i=1;
	for (i=1; i<=n; i++)
	{
		acc += mangoldt_lambda (i);
	}

	return acc;
}

int main ()
{
	int i;

	int nmax = 1000;
	
	long double sum = 0.0;
	int scale = 1;
	for (i=1; i<nmax; i++)
	{
		// sum = chebyshev (i);
		sum += mangoldt_lambda(i);

		if (0 == i%scale) 
		{
			printf ("%d	%d\n", i, mert);
			fflush (stdout);
			// if (i > scale*300) scale *= 2;
		}
	}
}

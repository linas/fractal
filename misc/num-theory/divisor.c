
/*
 * divisor.c
 *
 * Graphs of series involving divisor function.
 * These are ordinary x-y plots; output is 
 * ascii list of x-y values
 *
 * Linas Vepstas July 2006
 */

#include <math.h>
#include <stdio.h>

#include "modular.h"

long double divisor_series (long double x)
{
	long double acc = 0.0;

	long double xp = 1.0;
	int n=1;
	while (1)
	{
		long double term = xp * divisor (n);
		acc += term;

		if (term < 1.0e-20*acc) break;
		xp *= x;
		n++;
	}

	return acc;
}

int main ()
{
	int i;

	int nmax = 25000;
	int scale = 10;

	int sum = 0;
	for (i=1; i<=nmax; i++)
	// for (i=nmax; i>0; i--)
	{
		long double x = ((double) i)/((double) nmax);

		// long double y = divisor_series (x);

		// printf ("%d	%Lg	%26.18Lg\n", i, x, y);
		
		int d = divisor (i);
		sum += d;
		if (0 == i%scale)
		{
			long double lead = sum;
			lead -= i * logl(i) - i + 0.1544L*i;
			// lead -= i;
			// int d = sigma (i,1);
			printf ("%d	%d	%26.18Lg\n", i, d, lead);
		}
	}
}

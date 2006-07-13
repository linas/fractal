
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

#include "harmonic.h"
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

	int nmax = 1000000000;
	int scale = 1000000;

	int sum = 0;
	long double inf = 1.0e100;
	long double sup = -1.0e100;
	for (i=1; i<=nmax; i++)
	// for (i=nmax; i>0; i--)
	{
		long double x = ((double) i)/((double) nmax);

		// long double y = divisor_series (x);

		// printf ("%d	%Lg	%26.18Lg\n", i, x, y);
		
		int d = divisor (i);
		sum += d;
		long double del = sum - (i * logl(i) + i*(2.0*M_GAMMA -1.0L));
		if (sup < del) sup=del;
		if (inf > del) inf=del;
		if (0 == i%scale|| i==1)
		{
			long double lead = sum;
			lead = i * logl(i) + 2.0*i*(M_GAMMA -1.0L);
			// lead -= i;
			// int d = sigma (i,1);
			// printf ("%d	%d	%26.18Lg\n", i, d, sum-lead);
			printf ("%d	%d	%26.18Lg\n", i, d, inf);
			printf ("%d	%d	%26.18Lg\n", i, d, sup);
			inf= 1.0e100;
			sup=-1.0e100;
			fflush (stdout);
		}
	}
}

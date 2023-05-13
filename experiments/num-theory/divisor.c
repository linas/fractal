
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
#include "moebius.h"
#include "modular.h"

int main ()
{
	long long int i;

	long long int nmax = 10000000000ll;
	int scale = 10000000;

	nmax = 250;
	scale=1;

	long double sum = 0.0L;
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
			d = sigma (i,3);
			printf ("%Ld	%d	%26.18Lg\n", i, d, sum);
			// printf ("%Ld	%d	%26.18Lg\n", i, d, sum-lead);
			// printf ("%Ld	%d	%26.18Lg\n", i, d, inf);
			// printf ("%Ld	%d	%26.18Lg\n", i, d, sup);
			inf= 1.0e100;
			sup=-1.0e100;
			fflush (stdout);
		}
	}
}

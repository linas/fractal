
/*
 * ndist.c
 *
 * graph of totient of dyadics on unit interval
 *
 * Linas Vepstas December 2004
 */

#include <stdio.h>
#include "gcf.h"
#include "totient.h"

int main ()
{
	int i;

	int nmax = 720;

	int sum = 0;
	for (i=1; i<nmax; i++)
	{
		double x = ((double) i)/((double) nmax);

		int g = gcf32 (i, nmax);
		int n = i / g;
		int y = totient_phi (n);

		sum += y;

		printf ("%d	%d	%g	%d	%d\n", i, n, x, y, sum);
	}
}

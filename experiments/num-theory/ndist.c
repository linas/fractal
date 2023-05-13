
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

	// int nmax = 720;
	// int nmax = 719;
	int nmax = 512;

	double sum = 0.0;
	for (i=1; i<nmax; i++)
	{
		double x = ((double) i)/((double) nmax);

		int g = gcf32 (i, nmax);
		int nume = i / g;
		int deno = nmax / g;
		double y = totient_phi (nume);
		y /= (double) nume;

		sum += y/((double) nmax);

		printf ("%d	%d	%g	%g	%g\n", i, nume, x, y, sum);
	}
}

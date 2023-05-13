/*
 * pegs.C
 *
 * Measure the distance to the nearest peg in a grid.
 *
 * Linas Vepstas March 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "gcf.h"

int
main (int argc, char * argv[]) 
{
	int i,p,q;

	int nmax = atoi (argv[1]);

	for (i=1; i<nmax; i++)
	{
		q = nmax;
		int cf = gcf32(i,nmax);
		p = i/cf;
		q = nmax / cf;

		double x = ((double)p) / ((double) q);
		double dist = p*p+q*q;
		dist = sqrt(dist);
		printf ("%d\t%d\t%d\t%g\t%g\n", i, p, q, x, dist);
	}
}

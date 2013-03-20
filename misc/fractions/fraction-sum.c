/*
 * Sums of fractions, accidentally discovered in lignuistic data.
 *
 * Linas March 2013
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double * count(int nbins, double binsz, int maxcnt)
{
	double * bins = (double*) malloc (nbins*sizeof(double));

	for (int i = 1; i< maxcnt; i++)
	{
		double deno = 1.0 / ((double) i);
		for (int j=1; j<maxcnt; j++)
		{
			double x = j * deno;

			int ib = (int) floor(x / binsz);
			if (nbins <= ib) break;
			bins[ib] += deno;
		}
	}

	return bins;
}

int 
main (int argc, char * argv[])
{
	double binsz = 0.01;
	double * bins = count(10000, binsz, 1000);

	for (int i=0; i< 1000; i++)
	{
		printf("%d	%g	%g\n", i, i*binsz, bins[i]);
	}
}

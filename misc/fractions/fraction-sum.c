/*
 * Sums of fractions, accidentally discovered in lignuistic data.
 *
 * Linas March 2013
 */

#include <stdlib.h>
#include <math.h>

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
			bins[ib] += 1.0;
		}
	}

	return bins;
}

main ()
{
}

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
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s <binsz> <maxcnt>\n", argv[0]);
		exit(1);
	}

	double binsz = 0.01;
	binsz = atof(argv[1]);
	int maxcnt = atoi(argv[2]);
	double * bins = count(10000, binsz, maxcnt);

	printf("#\n# Data generated with binsz = %g maxcnt=%d\n#\n", binsz, maxcnt);

	double accum = 0.0;
	for (int i=0; i< 10000; i++)
	{
		// accum += bins[i] - 0.4;
		// accum += bins[i] - 0.43;
		accum += bins[i];
		printf("%d	%g	%g	%g\n", i, i*binsz, bins[i], accum);
	}
}

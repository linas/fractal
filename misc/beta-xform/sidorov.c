/*
 * sidorov.c
 *
 * Understand bit-sequence mappings.  The expander and compresor
 * functions. This time with alternative expansions, per sidorov.
 *
 * Linas Vepstas Dec 2017; Sept 2020
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Note that beta = 2*K
double sdr(double y, double K)
{
	// Iterate on y using mashed Bernoulli, and extract symbol dynamics
	char nbits[50];
	for (int i=0; i<50; i++)
	{
		if (0.5 <= y)
		{
			y -= 0.5;
			nbits[i] = 1;
		}
		else nbits[i] = 0;
		y *= 2.0*K;
	}

	// Reconstruct x in a ordinary bernoulli sequence.
	double acc = 1.0e-30;
	for (int i=0; i<50; i++)
	{
		acc *= 0.5;
		if (nbits[50-i-1])
		{
			acc += 0.5;
		}
	}

	return acc;
}

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit (1);
	}
	double lam = atof(argv[1]);

	int npts = 163;
	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		double z = sdr (x, lam);
		printf("%d	%g	%g\n", i, x, z);
	}
}

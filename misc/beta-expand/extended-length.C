/*
 * extended-length.C
 * Compute the extended length via the density histogram.
 *
 * Linas Vepstas Oct 2020
 */

#include <math.h>
#include "extended-density.C"

int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K nbits\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);
	int nbits = atoi(argv[2]);
	// int maxdepth = atoi(argv[3]);

#define NBINS 1000
	double histo[NBINS];
	extended_measure(2.0*Kay, histo, NBINS, nbits);

	int em = emrun(Kay);
	double beta = 2.0 * Kay;
	double betam = pow(beta, em);
	double lower = 1.0 / (2.0*beta);
	double upper = lower * (1.0 + 1.0/betam);

#define PHI (0.5 * (sqrt(5.0) + 1.0))
#define SCALE (0.5 * (PHI + 1.0))

	double prob = 0.0;
	for (int i=0; i<NBINS; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NBINS);
		x *= SCALE;

		if (lower < x and x <= upper) prob += histo[i];
		if (lower < x-0.5 and x-0.5 <= upper) prob += histo[i];
		if (lower < x-1.0 and x-1.0 <= upper) prob += histo[i];
	}
	prob /= NBINS;
	printf("duude beta=%g prob=%g len=%g\n", beta, prob, 1.0/prob);
}

// ================================================================

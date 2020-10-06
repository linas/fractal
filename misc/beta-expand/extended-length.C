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
		fprintf(stderr, "Usage: %s maxdepth nsamples\n", argv[0]);
		exit (1);
	}
	int nbits = 95;
	int maxdepth = atoi(argv[1]);
	int nsamples = atoi(argv[2]);

#define NBINS 5000
	printf("#\n# Extended slice\n");
	printf("#\n# nbits=%d nbins=%d\n", nbits, NBINS);
	printf("# maxdepth=%d nsamples=%d\n#\n", maxdepth, nsamples);
	printf("# Columns:\n# j beta prob length fraction\n");
	fflush(stdout);

#define NPTS 231
	for (int j=0; j<NPTS; j++)
	{
		double y = (((double) j) + 0.5)/ ((double) NPTS);
		double beta = 1.0 + y;
		// beta = 1.5 - 0.5*y;
		beta = 2.0 - 0.5*y;
		double Kay = 0.5*beta;

		double histo[NBINS];
		extended_measure(beta, maxdepth, nsamples, histo, NBINS, nbits);

		int em = emrun(Kay);
		double betam = pow(beta, em);
		double lower = 1.0 / (2.0*beta);
		double upper = lower * (1.0 + 1.0/betam);

#define PHI (0.5 * (sqrt(5.0) + 1.0))
#define SCALE (0.5 * (PHI + 1.0))

		double prob = 0.0;
		int bcnt = 0;
		for (int i=0; i<NBINS; i++)
		{
			double x = (((double) i) + 0.5)/ ((double) NBINS);
			x *= SCALE;

			if (lower < x and x <= upper) { prob += histo[i]; bcnt++; }
			if (lower < x-0.5 and x-0.5 <= upper) { prob += histo[i]; bcnt++; }
			if (lower < x-1.0 and x-1.0 <= upper) { prob += histo[i]; bcnt++; }
		}
		prob /= NBINS;
		double frac = ((double) bcnt) / NBINS;
		printf("%d	%g	%g	%g	%g\n", j, beta, prob, 1.0/prob, frac);
		fflush(stdout);
	}
}

// ================================================================

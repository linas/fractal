/*
 * extended-measure.C
 * Compute the extended measure via histogram.
 *
 * Linas Vepstas Oct 2020
 */

#include <math.h>
#include "extended-density.C"

int main (int argc, char* argv[])
{
	if (argc < 5)
	{
		fprintf(stderr, "Usage: %s K nbits maxdepth nsamples\n", argv[0]);
		exit (1);
	}
	double Kay = atof(argv[1]);
	int nbits = atoi(argv[2]);
	int maxdepth = atoi(argv[3]);
	int nsamples = atoi(argv[4]);

#define NBINS 400

	printf("#\n# Extended measure\n");
	printf("#\n# K=%g beta=%g depth=%d nsamp=%d bits=%d\n#\n", Kay, 2.0*Kay,
	       maxdepth, nsamples, nbits);
	fflush(stdout);

	double histo[NBINS];
	extended_measure(2.0*Kay, maxdepth, nsamples, histo, NBINS, nbits);

	double histbase[NBINS];
	parry_measure(2.0*Kay, nsamples*(1<<maxdepth), histbase, NBINS, nbits);

#define PHI (0.5 * (sqrt(5.0) + 1.0))
#define SCALE (0.5 * (PHI + 1.0))

	for (int i=0; i<NBINS; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NBINS);
		double y = x * SCALE;
		printf("%d	%g	%g	%g	%g\n", i, y, histo[i], x, histbase[i]);
	}
}

// ================================================================

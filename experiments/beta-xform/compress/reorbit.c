/*
 * reorbit.c
 *
 * Compressor, re-coded. Histogram.
 *
 * Linas Vepstas Jan 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "compress.c"

int main (int argc, char* argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "Usage: %s beta npts\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int npts = atoi(argv[2]);

#define NHISTO 1000
	double histo[NHISTO];
	for (int i=0; i< NHISTO; i++) histo[i] = 0.0;

	double scale = 1.0;

	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		double y = cpr(x, 0.5*beta);

		for (int j=0; j<NBITS-12; j++)
		{
			int n = floor(NHISTO * y / scale);
			if (n < 0) n=0;
			if (NHISTO <= n) n = NHISTO-1;
			histo[n] += ((double) NHISTO) / ((double) npts);

			y *= 2.0;
			if (1.0 < y) y -= 1.0;
		}
	}

	printf("#\n# Re-orbit: beta=%g  npts=%d\n#\n", beta, npts);
	double sum = 0.0;
	for (int i=0; i< NHISTO; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NHISTO);
		sum += histo[i] / ((double) NHISTO * (NBITS-12));
		printf("%d	%g	%g	%g\n", i, x, histo[i], sum);
	}

	fflush (stdout);
	return 0;
}

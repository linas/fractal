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
		fprintf(stderr, "Usage: %s from-beta to-beta\n", argv[0]);
		exit (1);
	}
	double from_beta = atof(argv[1]);
	double to_beta = atof(argv[2]);
	int npts = 1231234; // atoi(argv[2]);

#define NHISTO 1000
	double histo[NHISTO];
	for (int i=0; i< NHISTO; i++) histo[i] = 0.0;

	double scale = 1.0;

	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		double y = cpr(x, 0.5*from_beta);
		double z = pdr(y, 0.5*to_beta);

		int n = floor(NHISTO * z / scale);
		if (n < 0) n=0;
		if (NHISTO <= n) n = NHISTO-1;
		histo[n] += ((double) NHISTO) / ((double) npts);

#if 0
		for (int j=0; j<NBITS-12; j++)
		{
			int n = floor(NHISTO * y / scale);
			if (n < 0) n=0;
			if (NHISTO <= n) n = NHISTO-1;
			histo[n] += ((double) NHISTO) / ((double) npts);

			y *= 2.0;
			if (1.0 < y) y -= 1.0;
		}
#endif
	}

	printf("#\n# Re-orbit: from-beta=%g  to-beta=%g npts=%d\n#\n",
		from_beta, to_beta, npts);
	double sum = 0.0;
	for (int i=0; i< NHISTO; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NHISTO);
		// sum += histo[i] / ((double) NHISTO * (NBITS-12));
		sum += histo[i] / ((double) NHISTO);
		printf("%d	%g	%g	%g\n", i, x, histo[i], sum);
	}

	fflush (stdout);
	return 0;
}

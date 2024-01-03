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
	if (argc != 4)
	{
		fprintf(stderr, "Usage: %s from-beta to-beta npts\n", argv[0]);
		exit (1);
	}
	double from_beta = atof(argv[1]);
	double to_beta = atof(argv[2]);
	int npts = atoi(argv[3]);

#define NHISTO 1000
	double histo[NHISTO];
	for (int i=0; i< NHISTO; i++) histo[i] = 0.0;

	double scale = 1.0;

	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		double y = cpr(x, 0.5*from_beta);

#ifdef ADJOINT
		double z = pdr(y, 0.5*to_beta);

		int n = floor(NHISTO * z / scale);
		if (n < 0) n=0;
		if (NHISTO <= n) n = NHISTO-1;
		histo[n] += ((double) NHISTO) / ((double) npts);
#endif

#if 1
		y -= floor(y);
		if (0.5 < y) y -= 0.5;
		for (int j=0; j<NBITS-12; j++)
		{
			int n = floor(NHISTO * 2.0* y / scale);
			if (n < 0) n=0;
			if (NHISTO <= n) n = NHISTO-1;
			histo[n] += ((double) NHISTO) / ((double) npts * (NBITS-12));

			y *= to_beta;
			if (0.5 < y) y -= 0.5;
		}
#endif
	}

	printf("#\n# Re-orbit: from-beta=%g  to-beta=%g npts=%d\n#\n",
		from_beta, to_beta, npts);
	double sum = 0.0;
	for (int i=0; i< NHISTO; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NHISTO);
		sum += histo[i] / ((double) NHISTO);
		printf("%d	%g	%g	%g\n", i, x, histo[i], sum);
	}

	fflush (stdout);
	return 0;
}

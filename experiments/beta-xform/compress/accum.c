/*
 * accum.c
 *
 * Where does the expander function accumulate? A histogram.
 * I thought I did this once before, but I can't find the result.
 * I even recall looking at the figures...
 * Found it. It's in "../hessenberg/sampler.c"
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

	// double scale = pdr(1.0, 0.5*beta);
	double scale = cpr(1.0, 0.5*beta);

	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		// double y = pdr(x, 0.5*beta);
		double y = cpr(x, 0.5*beta);

		int n = floor(NHISTO * y / scale);
		if (n < 0) n=0;
		if (NHISTO <= n) n = NHISTO-1;
		histo[n] += ((double) NHISTO) / ((double) npts);
	}

	printf("#\n# Accum histogram beta=%g  npts=%d\n#\n", beta, npts);
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

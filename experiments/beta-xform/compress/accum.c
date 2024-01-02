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

// The expander. A binary bit sequence for x is obtained, and
// resummed using K.  Recall K=beta/2.
double pdr(double x, double Kay)
{
#define NBITS 50
	// Decompose x into a bit sequence.
	char nbits[NBITS];
	for (int i=0; i<NBITS; i++)
	{
		if (0.5 <= x)
		{
			x -= 0.5;
			nbits[i] = 1;
		}
		else nbits[i] = 0;
		x *= 2.0;
	}

	// Reconstruct x in a mashed bernoulli sequence.
	double acc = 1.0e-30;
	for (int i=0; i<NBITS; i++)
	{
		acc *= 1.0 / (2.0*Kay);
		if (nbits[NBITS-i-1])
			acc += 0.5;
	}

	return acc;
}

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s beta npts\n", argv[0]);
		exit (1);
	}
	double beta = atof(argv[1]);
	int npts = atoi(argv[2]);

#define NHISTO 1000
	double histo[NHISTO];
	for (int i=0; i< NHISTO; i++) histo[i] = 0.0;

	double scale = pdr(1.0, 0.5*beta);

	for (int i=0; i<npts; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) npts);
		double y = pdr(x, 0.5*beta);

		int n = floor(NHISTO * y / scale);
		if (n < 0) n=0;
		if (NHISTO <= n) n = NHISTO-1;
		histo[n] += ((double) NHISTO) / ((double) npts);
	}

	for (int i=0; i< NHISTO; i++)
	{
		double x = (((double) i) + 0.5)/ ((double) NHISTO);
		printf("%d	%g	%g\n", i, x, histo[i]);
	}

	fflush (stdout);
	return 0;
}

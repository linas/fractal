/*
 * witt.c
 * Witt-vector-like games
 * Currently failing. Fooo.
 *
 * March 2018
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double witt_xform(double x, double beta, int n)
{
	double prod = pow (beta, n);
	return x - floor(x*prod)/prod;
}

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit(1);
	}

	double K = atof(argv[1]);
	double beta = 2.0*K;

#define NPTS 1603
	double bins[NPTS];
	for (int i=0; i<NPTS; i++)
	{
		bins[i] = 0.0;
	}

#define NREPS 50
	double bn = beta;
	for (int n=0; n<NREPS; n++)
	{
		for (int i=0; i<NPTS; i++)
		{
			double x = ((double) i + 0.5) / ((double) NPTS);
			double y = bn*x - floor(x*bn);
			// double z = bn*witt_xform(x, beta, n+1);

			// printf("%d %d x=%g bn=%g y=%g\n", n, i, x, bn, y);
			int k = floor(y*NPTS);
			if (NPTS <= k) { k=NPTS-1; printf("ohhh no\n"); }
			bins[k]++;
		}
		bn *= beta;
	}

	double norm = 1.0 / ((double) NREPS);

	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double y = norm * bins[i];
		printf("%d	%g	%g\n", i, x, y);
	}
}

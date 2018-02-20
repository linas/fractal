/*
 * expan.c
 *
 * Numerical integration to obtain measure coefficients
 * Dec 2017
 * Feb 2018
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NOMAIN
#include "psi.c"

// Compute eigenfunction, recursively.
double reig(double x, double K, int niter)
{
	if (K < x) return 0.0;
	if (niter < 0)
	{
		// Approximate by a constant.
		return 1.0 / K;
	}
	double tkay = 2.0*K;

	// Short-cut. The other branch has vanishing contribution,
	if (K*(tkay-1.0) < x)
	{
		return reig(x/tkay, K, niter-1) / (lambda * tkay);
	}
	double sum = reig(x/tkay, K, niter-1);
	sum += reig(0.5 + x/tkay, K, niter-1);
	sum /= lambda * tkay;
	return sum;
}


/* Returns the n'th coefficient in the expansion of the invariant
 * measure. Computed by brute-force integration.
 */
double alpha_n(double K, int n, int npts)
{
	/* Get the lower, middle and upper bounds */
	n++; /* Off-by-one! */
	double lower = midpoints[lower_sequence[n]];
	double upper = midpoints[upper_sequence[n]];
	double middle = midpoints[n];

	/* integrate the lower arm */
	double delta = (middle-lower) / ((double) npts);
	double sumlo = 0.0;
	for (int i=0; i< npts; i++)
	{
		double x = lower + ((double) i + 0.5) * delta;
		sumlo +=
	}
}


int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s K\n", argv[0]);
		exit(1);
	}

	double K = atof(argv[1]);

	find_midpoints(K, MAXN);
	sequence_midpoints(K, MAXN);

	int npts = 10;
	double x = K;
	for (int i=0; i<npts; i++)
	{
		int bit = 0;
		if (0.5 < x) bit = 1;
		printf ("%d %g	%d	%g\n", i, x, bit, midpoints[i]);
		x = bmap(K, x);
	}
}

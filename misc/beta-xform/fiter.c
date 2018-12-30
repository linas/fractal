/*
 * fiter.c
 *
 * Attempt to fit the invariant measure, manually.
 * Feb 2018
 */

#include <stdio.h>
#include <stdlib.h>

#define NOMAIN
#include "psi.c"

/* Return the n'th wave function at the value of x. */
/* Differs from psi_n by normalization */
double chi_n(double x, double K, int n)
{
	// printf("psi ask for %d x=%g K=%g\n", n, x, K);
	if (K < x) return 0.0;
	if (0 == n)
	{
		return 1.0 / K;
		return 1.0 / sqrt(K);
	}

	/* Get the lower, middle and upper bounds */
	n++; /* Off-by-one! */
	double lower = midpoints[lower_sequence[n]];
	if (x < lower) return 0.0;
	double upper = midpoints[upper_sequence[n]];
	if (upper < x) return 0.0;

	double middle = midpoints[n];
	double norm = 1.0 / (middle - lower);
	norm += 1.0 / (upper - middle);
	norm = 1.0 / norm;
	// norm = sqrt(norm);
	if (x < middle)
		return norm / (middle - lower);

	return -norm / (upper - middle);
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

#if 0
	int npts = 10;
	double x = K;
	for (int i=0; i<npts; i++)
	{
		int bit = 0;
		if (0.5 < x) bit = 1;
		printf ("%d %g	%d	%g\n", i, x, bit, midpoints[i]);
		x = bmap(K, x);
	}
#endif
	int npts = 1601;
	for (int i=0; i<npts; i++)
	{
		double x = ((double) i + 0.5) / ((double) npts);
		printf("%d	%g", i, x);
		for (int j = 0; j<10; j++)
		{
			// double y = psi_n(x, K, j);
			double y = chi_n(x, K, j);
			printf("	%g", y);
		}
		printf("\n");
	}
}

/*
 * psimult.c
 * Multiplicative product of midpoints.
 *
 * February 2018
 */
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define NOMAIN
#include "psi.c"
#include "psibig.c"
#undef NOMAIN

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s maxn\n", argv[0]);
		exit(1);
	}
	int maxn = atoi(argv[1]);

	printf("#\n# Produce of %d midpoints\n#\n", maxn);

#define NPTS 1701
	double accum = 0.0;
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double K = 0.5 + 0.5*x;

		// find_midpoints(K);
		K = big_midpoints(K, 400, midpoints, maxn);
		sequence_midpoints(K, maxn);

		// product of midpoints
		double prod = 1.0;
		for (int m=1; m< maxn; m++)
		{
			prod *= 2.0 * midpoints[m] / K;
			// prod /= K;
		}
		accum += prod;
		printf("%d	%g	%g	%g\n", i, K, prod, accum);
	}
}

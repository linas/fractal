/*
 * walker.c
 * Graph walker
 *
 * Feb 2024
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "selfie.c"

// Walker. As currently written, this is correct only for
// beta=golden mean.
double walk(double beta, double eps, double x)
{
	double sum = 0.0;
	// int steps = 0;
	double pt = 0.5*beta;
	while (eps < fabs(x))
	{
		pt /= beta;
		if (0.5 < x)
		{
			// This loop is correct only for golden mean.
			// General case requires post-raise fall-downs.
			while (0.5*beta < pt + 0.5) pt /= beta;
			pt += 0.5;

			// Fall down once.
			pt /= beta;

			// Relocate
			sum += pt;
		}
		x *= 2.0;
		x -= floor(x);
	}
	return sum;
}

// -------------------------------------------------------

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: %s <index>\n", argv[0]);
		exit(1);
	}
	long aidx = atol(argv[1]);
	long idx = aidx;

	while (false == valid_gold_index(idx))
		idx++;

	if (idx != aidx)
		fprintf(stderr, "Warning: bump %ld to %ld\n", aidx, idx);

	double beta = golden_beta(idx);
	printf("# Beta = %f\n", beta);

#define NPTS 512
	double eps = 0.5 / ((double) NPTS);
	for (int i=0; i< NPTS; i++)
	{
		double x = ((double) i) / ((double) NPTS);
		double y = walk(beta, eps, x);

		printf("%d	%f	%g\n", i, x, y);
	}
}

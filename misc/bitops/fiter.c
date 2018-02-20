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

double bmap(double K, double x)
{
	if (0.5 < x)
		x -= 0.5;
	return 2.0*K*x;
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
			double y = psi_n(x, K, j);
			printf("	%g", y);
		}
		printf("\n");
	}
}

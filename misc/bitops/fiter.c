/*
 * fiter.c
 *
 * Attempt to fit the invariant measure, manually.
 * Feb 2018
 */

#include <stdio.h>
#include <stdlib.h>

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
	double x = K;

	int npts = 10;
	int bits[npts];

	for (int i=0; i<npts; i++)
	{
		int bit = 0;
		if (0.5 < x) bit = 1;
		printf ("%d %g	%d\n", i, x, bit);
		x = bmap(K, x);
	}
}

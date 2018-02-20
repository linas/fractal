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

	for (int i=0; i< 10; i++)
	{
		printf ("%d %g\n", i, x);
		x = bmap(K, x);
	}
}

/*
 * orbit.c
 * Graph of midpoint orbit
 *
 * February 2018
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NOMAIN
#include "psi.c"
#include "psibig.c"
#undef NOMAIN

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K npts\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	int npts = atoi(argv[2]);

	// find_midpoints(K);
	big_midpoints(K, 400, midpoints, npts);
	for (int i=0; i<npts; i++)
	{
		printf("%d	%g\n", i, midpoints[i]);
	}
}

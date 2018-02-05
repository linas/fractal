/*
 * psimom.c
 * Experiments with modelling the moments Moment matrix with actual
 * momments.
 *
 * February 2018
 */
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#define NOMAIN
#include "psileft.c"
#undef NOMAIN


int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K maxn\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	int maxn = atoi(argv[2]);

	printf("#\n# K=%g\n#\n", K);

	// find_midpoints(K);
	big_midpoints(K, 400, midpoints, MAXN);
	sequence_midpoints(K, MAXN);

#ifdef COLUMN_RATIOS
	// Well, its symmetric, so it doesn't matter: rows or cols.
	int n = 1;
	double prev = herm(K, n, 0);
	for (int m=1; m<maxn; m++)
	{
		double sym = herm(K, n, m);
		double rat = sym/prev;
		printf("%d	%d %g %g\n", n, m, sym, rat);
		fflush(stdout);
		prev = sym;
	}
#endif

#define HI_ACC_RATIOS
#ifdef HI_ACC_RATIOS
	int n = maxn;
	int m = maxn;
	// m = maxn*0.723;
	// m = maxn*0.411;
	// m = maxn*1.711;
	double prev = herm(K, n, m-1);
	double sym = herm(K, n, m);
	double rat = sym/prev;
	printf("%g	%d	%d %g %g\n", K, n, m, sym, rat);
#endif

}

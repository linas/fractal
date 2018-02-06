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

#define NOMOMAIN
#include "psileft.c"
#undef NOMOMAIN


int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K maxn\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	int maxn = atoi(argv[2]);

K = 0.81 + K*0.001 * 0.19;

	printf("#\n# K=%g\n#\n", K);

	// find_midpoints(K);
	big_midpoints(K, 400, midpoints, MAXN);
	sequence_midpoints(K, MAXN);

	// Get a high-accuracy ratio
	int n = maxn;
	int m = maxn;
	double prev = herm(K, n, m-1);
	double sym = herm(K, n, m);
	double rat = sym/prev;

	double pnt = pow(rat, n+m+1);
	double lng = sym/pnt;
	// printf("# %g	%d	%d sym=%g rat=%g f=%g\n#\n", K, n, m, sym, rat, lng);
	printf("%g	%d	%d %g %g %g\n", K, n, m, sym, rat, lng);


// #define PRINT_MATRIX
#ifdef PRINT_MATRIX
	for (int n=0; n<maxn; n++)
	{
		for (int m=0; m<maxn; m++)
		{
			double her = herm(K, n, m);
			double pnt = pow(rat, n+m+1);
			// double asy = her/ (pnt * lng);
			double rem = her - pnt*lng;
			// printf("[%d	%d] =	%g	%g	%g	%g\n", n, m, her, pnt, asy, rem);
			printf("[%d	%d] =	%g	%g\n", n, m, her, rem);
		}
		printf(" ---------------\n");
	}
#endif

}

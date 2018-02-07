/*
 * psivec.c
 * Experiments with vectors pushed through the polynomial matrixes.
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
	if (argc < 4)
	{
		fprintf(stderr, "Usage: %s K w maxn\n", argv[0]);
		exit(1);
	}
	double K = atof(argv[1]);
	double w = atof(argv[2]);
	int maxn = atoi(argv[3]);

	printf("#\n# K=%g	w=%g\n#\n", K, w);

	// find_midpoints(K);
	big_midpoints(K, 400, midpoints, MAXN);
	sequence_midpoints(K, MAXN);

	// create vector
	double wvec[maxn];
	double wn = 1.0;
	for (int n=0; n<maxn; n++)
	{
		wvec[n] = wn;
		wn *= w;
	}

	// multiply by p inverse transpose
	double pitw[maxn];
	for (int n=0; n<maxn; n++)
	{
		double acc = 0.0;
		for (int m=0; m<maxn; m++)
		{
			double pit = rinverse(K, m, n); // inverse transpose
			acc += pit * wvec[m];
		}
		pitw[n] = acc;
	}

	// renormalize
	double ren = pitw[0];
	for (int n=0; n<maxn; n++)
	{
		pitw[n] /= ren;
	}

	// multiply by hess
	double ha[maxn];
	for (int n=0; n<maxn; n++)
	{
		double acc = 0.0;
		for (int m=0; m<maxn; m++)
		{
			double he = hess(K, n, m);
			acc += he * pitw[m];
		}
		ha[n] = acc;
	}

	for (int n=0; n<maxn; n++)
	{
		double rat = ha[n] / pitw[n];
		printf("# %d	%g	%g	%g	%g\n", n, wvec[n], pitw[n], ha[n], rat);
	}

#ifdef PRINT_INVARIANT
#define NPTS 801
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double y = 0.0; 
		for (int n=0; n<maxn; n++)
		{
			y += ha[n] * psi_n(x, K, n);
		}
		printf("%d	%g	%g\n", i, x, y);
	}
#endif

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

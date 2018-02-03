/*
 * psileft.c
 * Compute the Bergman (left-hand) polynomial coefficients.
 * The bergman.C file does similar, but only uses Z.
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


/** Evaluate the entries in the Bergman vector at some fixed z */
void bergman_vect(double Kay, complex z, complex* poly, int maxn)
{
	poly[0] = 1.0;
	for (int n=1; n<maxn; n++)
	{
		complex acc = z * poly[n-1];
		for (int k=0; k<n; k++)
		{
			acc -= hess(Kay, k, n-1) * poly[k];
		}
		acc /= hess(Kay, n, n-1);
		poly[n] = acc;
	}
}

/** Just like above, but without the z values */
double bergman_oper(double Kay, int n, int j)
{
	if (0 == n && 0 == j) return 1.0;
	if (n < j) return 0.0;
	if (j < 0) return 0.0;

	double acc = bergman_oper(Kay, n, j-1);
	
	return acc;
}

complex eval_poly(complex z, complex* poly, int maxn)
{
	complex zn = 1;
	complex acc = 0.0;
	for (int n=0; n<maxn; n++)
	{
		acc += zn * poly[n];
		zn *= z;
	}
	return acc;
}

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

	complex poly[MAXN];
	complex z = 1.0 / (2.0 * 0.8);
	z = 1.0;
	z = 0.5;
	z = 1.0 / (2.0 * K);
	bergman_vect(K, z, poly, MAXN);

	for (int i=0; i<maxn; i++)
	{
		printf("%d	%g	%g\n", i, creal(poly[i]), cimag(poly[i]));
	}
}

/*
 * reigen.c
 *
 * Recursve eigenfunctions for the undershift
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double reig(double x, double K, double lambda, int niter)
{
	if (niter < 0) return 1.0;
	if (K < x) return 0.0;
	if (K*(2.0*K-1.0) < x)
	{
		return reig(0.5*x/K, K, lambda, niter-1) / (2.0*lambda*K);
	}
	double sum = reig(0.5*x/K, K, lambda, niter-1);
	sum += reig(0.5 + 0.5*x/K, K, lambda, niter-1);
	sum /= (2.0*lambda*K);
	return sum;
}


int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K w\n", argv[0]);
		exit (1);
	}
	double K = atof(argv[1]);
	double lambda = atof(argv[2]);

	printf("#\n# K = %g lambda = %g\n#\n", K, lambda);

#define NPTS 1803
// #define NRECU 40
#define NRECU 20
	double acc = 0.0;
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i) / ((double) NPTS);
		double y = reig(x, K, lambda, NRECU);
		acc += y;
	}
	acc /= NPTS;

	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i) / ((double) NPTS);
		double y = reig(x, K, lambda, NRECU);
		y /= acc;
		printf("%d	%g %g\n", i, x, y);
	}
}

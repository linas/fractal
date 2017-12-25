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
	if (K < x) return 0.0;
	if (niter < 0) return 1.0;

	double tkay = 2.0*K;
	if (K*(tkay-1.0) < x)
	{
		return reig(x/tkay, K, lambda, niter-1) / (lambda * tkay);
	}
	double sum = reig(x/tkay, K, lambda, niter-1);
	sum += reig(0.5 + x/tkay, K, lambda, niter-1);
	sum /= lambda * tkay;
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

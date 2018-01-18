/*
 * complex.c
 *
 * Recursve complex-valued eigenfunctions for the undershift 
 * Jan 2018
 */

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NPTS 1803
double hits[NPTS];

// Compute eigenfunction, recursively.
double complex reig(double x, double K, double complex lambda, int niter)
{
	if (niter < 0)
	{
		unsigned int bin = floor (x*NPTS);
		hits[bin] += 1.0;

		if (K < x) return 0.0;

		// Approximate by something that integrates to zero.
		double re = 1.0;
		if (0.5*K < x) re = -1.0;
		double im = 1.0;
      if (0.25*K < x && x < 0.75*K) im = -1.0;

		return re + I*im;
	}
	if (K < x) return 0.0;

	double tkay = 2.0*K;

	// Short-cut. The other branch has vanishing contribution,
	if (K*(tkay-1.0) < x)
	{
		return reig(x/tkay, K, lambda, niter-1) / (lambda * tkay);
	}
	double complex sum = reig(x/tkay, K, lambda, niter-1);
	sum += reig(0.5 + x/tkay, K, lambda, niter-1);
	sum /= lambda * tkay;
	return sum;
}


int main (int argc, char* argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Usage: %s K t\n", argv[0]);
		exit (1);
	}
	double K = atof(argv[1]);
	double theta = atof(argv[2]);
	double complex lambda = cexp (I*2.0*M_PI*theta) / (2.0*K);

	printf("#\n# K = %g theta = %g lambda = %g + I %g \n#\n",
	       K, theta, creal(lambda), cimag(lambda));

	for (int i=0; i<NPTS; i++)
		hits[i] = 0.0;

// #define NRECU 30
#define NRECU 20

	// Compute an eigenfunction, recursively.
	for (int i=0; i<NPTS; i++)
	{
		double x = ((double) i + 0.5) / ((double) NPTS);
		double complex y = reig(x, K, lambda, NRECU);

		printf("%d	%g	%g	%g\n", i, x, creal(y), cimag(y));
	}
}

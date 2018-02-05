/*
 * moments.c
 *
 * Moments of the downshift.
 * Feb 2018
 */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Compute eigenfunction (invariant measure), recursively.
double reig(double x, double K, int niter)
{
	if (niter < 0)
	{
		return 1.0 / K;
	}
	if (K < x) return 0.0;

	double tkay = 2.0*K;

	// Short-cut. The other branch has vanishing contribution,
	if (K*(tkay-1.0) < x)
	{
		return reig(x/tkay, K,  niter-1) / tkay;
	}
	double sum = reig(x/tkay, K, niter-1);
	sum += reig(0.5 + x/tkay, K, niter-1);
	sum /= tkay;
	return sum;
}


int main (int argc, char* argv[])
{
	if (argc < 5)
	{
		fprintf(stderr, "Usage: %s K n-moments npts nrecur\n", argv[0]);
		exit (1);
	}
	double K = atof(argv[1]);
	int nmom = atoi(argv[2]);
	int npts = atoi(argv[3]);
	int nrecu = atoi(argv[4]);

	printf("#\n# K = %g npts = %d\n#\n", K, npts);

	
	// Compute invariant measure
	double mu[npts];
	double sum = 0.0;
	for (int i=0; i<npts; i++)
	{
		double x = ((double) i + 0.5) / ((double) npts);
		x *= K;
		mu[i] = reigen(x, K, nrecu);
		sum += mu[i]*mu[i];
	}

	// Renormalize
	sum /= npts;
	sum = 1.0 / sqrt(sum);
	for (int i=0; i<npts; i++)
	{
		mu[i] *= sum;
	}

	for (int k=0; k<nmom; k++)
	{
		double mom = 0.0;
		for (int i=0; i<npts; i++)
		{
			double x = ((double) i + 0.5) / ((double) npts);
			double xn = pow(x, k);
			mom += mu[i] * xn;
		}
		mom /= npts;
		printf("%d	%g\n", k, mom);
	}
}

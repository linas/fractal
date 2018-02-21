/*
 * expan.c
 *
 * Numerical integration to obtain measure coefficients.
 * This is a completely different, alternate take on psifp.c
 * and is much slower. It serves to cross-check results.
 *
 * Dec 2017
 * Feb 2018
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NOMAIN
#include "psi.c"
#include "psibig.c"

// Compute eigenfunction, recursively.
double reig(double x, double K, int niter)
{
	if (K < x) return 0.0;
	if (niter < 0)
	{
		// Approximate by a constant.
		return 1.0 / K;
	}
	double tkay = 2.0*K;

	// Short-cut. The other branch has vanishing contribution,
	if (K*(tkay-1.0) < x)
	{
		return reig(x/tkay, K, niter-1) / tkay;
	}
	double sum = reig(x/tkay, K, niter-1);
	sum += reig(0.5 + x/tkay, K, niter-1);
	sum /= tkay;
	return sum;
}


/* Returns the n'th coefficient in the expansion of the invariant
 * measure. Computed by brute-force integration.
 */
double alpha_n(double K, int n, int npts, int niter)
{
	if (0 == n)
	{
		/* integrate everything */
		double delta = K / ((double) npts);
		double sum = 0.0;
		for (int i=0; i< npts; i++)
		{
			double x = ((double) i + 0.5) * delta;
			sum += reig(x, K, niter);
		}
		sum *= delta;

		/* Normalize per the wavelet normalization */
		double norm = 1.0 / sqrt(K);
		sum *= norm;
		return sum;
	}

	/* Get the lower, middle and upper bounds */
	n++; /* Off-by-one! */
	double lower = midpoints[lower_sequence[n]];
	double upper = midpoints[upper_sequence[n]];
	double middle = midpoints[n];

	/* integrate the lower arm */
	double delta = (middle-lower) / ((double) npts);
	double sumlo = 0.0;
	for (int i=0; i< npts; i++)
	{
		double x = lower + ((double) i + 0.5) * delta;
		sumlo += reig(x, K, niter);
	}
	sumlo *= delta;

	/* integrate the upper arm */
	delta = (upper-middle) / ((double) npts);
	double sumhi = 0.0;
	for (int i=0; i< npts; i++)
	{
		double x = middle + ((double) i + 0.5) * delta;
		sumhi += reig(x, K, niter);
	}
	sumhi *= delta;

	/* Combine the arms, according to the wavelet size. */
	double sum = sumlo / (middle - lower);
	sum -= sumhi / (upper - middle);

	/* The wavelet normalization */
	double norm = 1.0 / (middle - lower);
	norm += 1.0 / (upper - middle);
	norm = 1.0 / sqrt(norm);

	sum *= norm;

	return sum;
}


int main(int argc, char* argv[])
{
	if (argc < 4)
	{
		fprintf(stderr, "Usage: %s K npts niter\n", argv[0]);
		exit(1);
	}

	double K = atof(argv[1]);
	int npts = atoi(argv[2]);
	int niter = atoi(argv[3]);

	printf("#\n# npts = %d depth=%d\n#\n", npts, niter);

	// Total integral is unity
	double norm = 0.0;
	for (int i = 0; i< npts; i++)
	{
		double x = ((double) i + 0.5) / ((double) npts);
		norm += reig(x, K, niter);
	}
	norm /= ((double) npts);
	printf("# tot integral = %g\n", norm);
	norm = 1.0 / norm;

	// find_midpoints(K, MAXN);
	big_midpoints(K, 400, midpoints, MAXN);
	sequence_midpoints(K, MAXN);

	for (int i=0; i<200; i++)
	{
		double cf = norm * alpha_n(K, i, npts, niter);
		if (cf < 2.0e-15) break;
		printf ("%d %g\n", i, cf);
		fflush(stdout);
	}
}

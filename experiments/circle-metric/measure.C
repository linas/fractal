/*
 * measure.C
 * Quick hack bin-count measure for selected points on the circlemap.
 *
 * History:
 * Linas Vepstas November 2023
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*-------------------------------------------------------------------*/
/*
 * Compute the measure for the circle map, for a selected location.
 * Print the normalized bincounts to stdout.
 *
 */

// Make sure that ITER_DEPTH is at least 4x larger than nbins,
// so that for chaotic orbits, each bin is hit maybe 4 times.
// This will avoid issues with the thresholding, below.
#define ISS_NBINS 1000
#define ISS_SETTLE_TIME 190
// #define ISS_ITER_DEPTH 1920       // Iteration depth
#define ISS_ITER_DEPTH 8111
// #define ISS_ITER_DEPTH 38111

void
circle_measure (double omega, double K, int itermax)
{
	// Initialize empty bins
	double bins[ISS_NBINS+1];
	for (int ib=0; ib<=ISS_NBINS; ib++)
		bins[ib] = 0;

	// Iterate, with random starts.
	for (int j=0; j<itermax/ISS_ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* First, we give a spin for 500 cycles, giving the non-chaotic
		 * parts a chance to phase-lock */
		for (int iter=0; iter<ISS_SETTLE_TIME; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
		}

		// OK, now, bin-count as we move along.
		// bin-counting is always modulo one, because that is
		// all that sine cares about.
		for (int iter=0; iter < ISS_ITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
			double yb = ISS_NBINS * (x - floor (x));
			int ib = (int) yb;
			bins[ib]++;
		}
	}

	// Get the total bincount.
	long totcnt = 0;
	for (int ib=0; ib<=ISS_NBINS; ib++)
		totcnt += bins[ib];

	// Dump to stdout
	printf("#\n# Cicle-map bincount, omega=%f K=%f\n", omega, K);
	printf("# nbins=%d totcnt=%ld\n", ISS_NBINS, totcnt);
	printf("# iter-depth=%d rand-start=%d itermax=%d\n#\n",
	       ISS_ITER_DEPTH, itermax/ISS_ITER_DEPTH, itermax);

	for (int ib=0; ib<=ISS_NBINS; ib++)
	{
		double x = ib + 0.5;
		x /= ISS_NBINS;
		double mu = bins[ib];
		mu /= totcnt;
		printf("%d	%f	%f\n", ib, x, mu);
	}
}

int main(int argc, char* argv[])
{
	circle_measure (0.1, 0.97, 10000);
}

/* --------------------------- END OF LIFE ------------------------- */

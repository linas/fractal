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

// #define ISS_SETTLE_TIME 1190

// Points below K=1/2pi have a vanishing Lyapunov const, and need a huge
// settle time. Laypunov const seems like its about 10e-4 ???
// .. Except this does not seem to matter?
// #define ISS_SETTLE_TIME 31190
#define ISS_SETTLE_TIME 190

#define LYA_DELTA 1e-5

void
circle_measure (double omega, double K, int nbins, int nstarts, int depth)
{
	// Initialize empty bins
	double bins[nbins+1];
	double lyap[nbins+1];
	for (int ib=0; ib<=nbins; ib++)
	{
		bins[ib] = 0;
		lyap[ib] = 0;
	}

	// Iterate, with random starts.
	for (int j=0; j<nstarts; j++)
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
		for (int iter=0; iter < depth; iter++)
		{
			double xnought = x;
			x += omega - K * sin (2.0 * M_PI * x);
			double yb = nbins * (x - floor (x));
			int ib = (int) yb;
			bins[ib]++;

			xnought += LYA_DELTA;
			xnought += omega - K * sin (2.0 * M_PI * xnought);
			double delta = fabs(x - xnought);
			double expo = log(delta / LYA_DELTA);
			lyap[ib] += expo;
		}
	}

	// Get the total bincount.
	long totcnt = 0;
	for (int ib=0; ib<nbins; ib++)
		totcnt += bins[ib];

	// Get average lyapunov for this bin
	for (int ib=0; ib<nbins; ib++)
		lyap[ib] /= bins[ib];

	// Dump to stdout
	printf("#\n# Cicle-map bincount, omega=%f K=%f\n", omega, K);
	printf("# nbins=%d totcnt=%ld\n", nbins, totcnt);
	printf("# iter-depth=%d rand-start=%d\n#\n", depth, nstarts);

	for (int ib=0; ib<nbins; ib++)
	{
		double x = ib + 0.5;
		x /= (double) nbins;
		double mu = bins[ib];
		mu /= (double) totcnt;
		mu *= (double) nbins;
		printf("%d	%f	%f	%f\n", ib, x, mu, lyap[ib]);
	}
}

int main(int argc, char* argv[])
{
	if (6 > argc)
	{
		fprintf(stderr,
			"Usage: %s <omega> <k> <nbins> <nstarts> <depth>\n", argv[0]);
		exit(1);
	}
	double omega = atof(argv[1]);
	double K = atof(argv[2]);
	double nbins = atoi(argv[3]);
	double starts = atoi(argv[4]);
	double depth = atoi(argv[5]);

	// circle_measure (0.1, 1.1, itermax);
	circle_measure (omega, K, nbins, starts, depth);
}

/* --------------------------- END OF LIFE ------------------------- */

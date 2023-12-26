/*
 * slice.c
 * Vertical slice through the measure
 *
 * December 2023
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*
 * The standard gelfond-parry measure.
 * (using our convention of 0 < x < beta/2)
 */
double gp_measure(double x, double beta)
{
#define NITER 180
	double norm = 0.0;
	double meas = 0.0;
	double midpoint = 0.5*beta;
	double obn = 1.0;
	for (int i=0; i<NITER; i++)
	{
		// Compare x to midpoint, and accumulate the measure.
		if (x < midpoint)
			meas += obn;

		// Accumulate the normalization
		norm += midpoint * obn;

		// Iterate the midpoint
		if (0.5 < midpoint)
			midpoint -= 0.5;
		midpoint *= beta;

		// Update 1/beta^n  for the next go-round
		obn /= beta;

		if (obn < 1e-15) break;
	}

	return meas/norm;
}

int main(int argc, char *argv[])
{
	double scale = atof(argv[1]);
#define NPTS 551
	double acc = 0.0;
	for (int i=0; i<NPTS; i++)
	{
		double v = (((double) i) + 0.5) / ((double) NPTS);
		double beta = 1.0 + v;
		double x = scale * 0.5*beta;

		double mu = gp_measure(x, beta);
		acc += mu / ((double) NPTS);

		printf("%d	%g	%g	%g\n", i, beta, mu, acc);
	}

#ifdef HORIZONTAL_SLICE
	// Conventional slice for fixed beta
	// This is just a sanity check of the code.
	double beta = atof(argv[1]);
	double acc = 0.0;
	for (int i=0; i<NPTS; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NPTS);

		double mu = gp_measure(x, beta);
		acc += mu / ((double) NPTS);

		printf("%d	%g	%g	%g\n", i, x, mu, acc);
	}
#endif
}

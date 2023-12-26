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
 * The standard gelfond-aprry measure.
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
		// First, compare x to midpoint, and accumulate the measure.
		if (x < midpoint)
			meas += obn;

		// Next, accumulate the normalization and iterate the midpoint
		if (0.5 < midpoint)
		{
			midpoint -= 0.5;
			norm += obn;
		}
		else
		midpoint *= beta;

		// Update 1/beta^n  for the next go-round
		obn /= beta;

		if (obn < 1e-15) break;
	}

	return meas/norm;
}

int main(int argc, char *argv[])
{
	double beta = atof(argv[1]);

	// Conventional slice for fixed beta
#define NPTS 551
	double acc = 0.0;
	for (int i=0; i<NPTS; i++)
	{
		double x = (((double) i) + 0.5) / ((double) NPTS);

		double mu = gp_measure(x, beta);
		acc += mu / ((double) NPTS);

		printf("%d	%g	%g	%g\n", i, x, mu, acc);
	}
}

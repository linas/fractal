/*
 * poincare.C
 *
 * Brand-new recreation of the circle-map stuff.
 * The original version is in generate/circle.C
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

double sign(double x)
{
	if (x < 0.0) return -1.0;
	if (0.0 < x) return 1.0;
	return 0.0;
}

// the downshift map, with half of an S-curve for the middle segment
// With the poincare convention K/beta will got left right, while
// epsilon will increase from bottom upwards.
double ess_island(double x, double K, double epsilon)
{
	epsilon *= 0.20;  // Epsilon runs from 0 to 0.2
	K = 0.5 + 0.5*K;  // K runs from 1/2 to 1

	double beta = 2.0 * K;
	if (0.5+epsilon <= x)
	{
		return beta * (x - 0.5);
	}
	if (x < 0.5-epsilon)
	{
		return beta*x;
	}

	// om runs from  -1 to 1;
	double om = (x - 0.5) / epsilon;

	// kos runs from -1 to 1
	// double kos = pow(om, 3);

	// With this sign convention, the map is continuous, and
	// has an S-curve shape.
	// double kos = om;
	// double kos = sign(om) * om*om;
	// double kos = om*om*om;
	// double kos = sign(om) * om*om*om*om;
	double kos = om*om*om*om*om;

	// With this sign convention, the middle segment is increasing,
	// and the map consists of three disjoint segments.
	// double kos = -om;
	// double kos = -sign(om) * om*om;
	// double kos = -om*om*om;
	// double kos = -sign(om) * om*om*om*om;
	// double kos = -om*om*om*om*om;

	return 0.25*beta - beta*(0.25-epsilon) * kos;
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the circle map
 */

// #define EPSILON 0.001
// #define EPSILON 0.009
#define EPSILON 0.009
#define SETTLE_TIME 21
#define RSAMP 1400

double
recurrance_time (double omeg, double Kba, int itermax,
                 double (func)(double, double, double) )

{
	double	x, y;
	double	xpoint;
	int		j, iter;
	long		num_recurs, time_recur=0;

	num_recurs = 0;
	for (j=0; j<itermax; j++)
	{
		double t = rand();
		t /= RAND_MAX;

// 800 x 800 pixels -- add some jitter.
// Except it needs to be correctly normalized, basd on magnification..!
// which we don't have available here.
// #define JITTER ((double) (800*2*2))
#define JITTER ((double) (800))
double omega = omeg + (t-0.5)/JITTER;
		t = rand();
		t /= RAND_MAX;
double Kbar = Kba + (t-0.5)/JITTER;
		t = rand();
		t /= RAND_MAX;
		x = t;

		// First, we give a spin for a few cycles, giving the
		// non-chaotic parts a chance to phase-lock.
		for (iter=0; iter<SETTLE_TIME; iter++)
		{
			x = func(x, omega, Kbar);
		}

		// OK, now, we begin to measure the average amount of time
		// to recur.
		xpoint = x;
		long ptime = 0;
		for (iter=0; iter < RSAMP; iter++)
		{
			x = func(x, omega, Kbar);
			y = fabs (x-xpoint);
			y -= floor (y);
			if (y < EPSILON)
			{
				num_recurs ++;
				ptime = iter;
			}
		}
		time_recur += ptime;
	}

	/* x is the (normalized) number of cycles to reach recurrance */
	x = (double) time_recur / ((double)num_recurs);

	return x;
}

static double hard_gram(double omega, double Kbar, int itermax, double param)
{
	return recurrance_time(omega, Kbar, itermax, ess_island);
}

DECL_MAKE_HEIGHT (hard_gram);

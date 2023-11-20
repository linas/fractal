/*
 * circle.C
 *
 * FUNCTION:
 * Circle map -- all new, all improved -- do it all over again.
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 * more stuff -- January 2000
 * more stuff -- October 2004
 * more stuff -- January 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/*-------------------------------------------------------------------*/
/*
 * This routine computes average winding number taken by
 * circle map iterator.
 */
static double winding_number (double omega, double K, int itermax)
{
	int cnt=0;
	double start=0.0, end=0.0;

#define ITER_DEPTH 150  // Number of iteration steps.
	for (int j=0; j<itermax/ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;
		start += x;

		/* OK, now start iterating the circle map */
		for (int iter=0; iter < ITER_DEPTH; iter++) {
			x += omega - K * sin (2.0 * M_PI * x);
			cnt ++;
		}
		end += x;
	}

	double wind = (end-start) / ((double) cnt);
	return wind;
}

/*-------------------------------------------------------------------*/
/*
 * This routine computes average winding number taken by
 * circle map iterator. -- subject to noise
 */
static double noisy_winding_number (double omega, double K, int itermax, double noise)
{
	double	x=0.0;
	int		iter;
	int cnt=0;

  	/* OK, now start iterating the circle map */
  	for (iter=0; iter < itermax; iter++) {
		x += omega - K * sin (2.0 * M_PI * x);
		cnt ++;

		/* white noise, equi-distributed, sharp cutoff */
		double t = rand();
		t /= RAND_MAX;
		// t -= 0.5;
		x += noise*t;
  	}

	x /= ((double) cnt);
	return x;
}

/*-------------------------------------------------------------------*/
/*
 * This routine computes average root-mean-square winding number
 * taken by circle map iterator.
 */
static double rms_winding_number (double omega, double K, int itermax)
{
	double	x=0.0, sq=0.0;
	int		iter,j;
	int cnt=0;
	double start=0.0, end=0.0;

	for (j=0; j<itermax/ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		x = t;
		start += x;

		/* OK, now start iterating the circle map */
		for (iter=0; iter < ITER_DEPTH; iter++) {
			x += omega - K * sin (2.0 * M_PI * x);
			sq += (x-t)*(x-t);
			t = x;
			cnt ++;
		}
		end += x;
	}

	// x = sqrt (sq) / (end-start);
	x = sqrt (sq) / ((double) cnt);
	// if (K != 0.0) x /= K;
	return x;
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the circle map
 */

#define EPSILON  	0.003
#define SETTLE_TIME 	90
#define RITER_DEPTH 500       // Iteration depth

double
circle_poincare_recurrance_time (double omega, double K, int itermax)

{
	double	x, y;
	double	xpoint;
	int		j, iter;
	long		num_recurs, time_recur=0;

  	num_recurs = 0;
	for (j=0; j<itermax/RITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		x = t;

		/* First, we give a spin for 500 cycles, giving the non-chaotic
		 * parts a chance to phase-lock */
		for (iter=0; iter<SETTLE_TIME; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
		}

		/* OK, now, we begin to measure the average amount of time to recur */
		/* (note that we don't have todo += with iter, since its already a running sum). */
		xpoint = x;
		long ptime = 0;
		for (iter=0; iter < RITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
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

/*-------------------------------------------------------------------*/
/*
 * Compute the Laplacian of the orbit (the charge) for the circle map.
 * This computes the series average for the four-point discrete
 * Laplacian for neighboring orbits.
 *
 * That is, consider an orbit (a time series) at fixed (K, omega). Then
 * consider the orbit for a nearby neighbor (K+deltaK, omega+deltaomega)
 * In mode-locked regions, we expect these to move together. At the
 * boundaries, we expect these to be uncorellated.
 *
 * At each point in the time-series, take the four-point Laplacian, viz:
 * 4x_n(K, o) - x_n(K-dK, o) - x_n(K+dK, o) - x_n(K, o-do) - x_n(K, o+do)
 * and then just average this together over all n.
 *
 * Except this is a stupid idea; neighboring pixels tend to have
 * opposite signs, when the are on opposite sides of an edge. Instead,
 * use the metric variants, below.
 */

// #define LAP_SETTLE_TIME 	90
#define LAP_SETTLE_TIME 	0

// Iteration depth
// #define LAP_ITER_DEPTH 50 // Ugh.
#define LAP_ITER_DEPTH 1000

double
circle_laplacian (double omega, double K, int itermax, double param)

{
	// Sample offsets. We want this to be the distance to the neighboring
	// pixel, more or less. User-specified. A hard-coded 0.001 is not a
	// bad place to start.
	double delta_K = param; // 0.001;
	double delta_omega = param; // 0.001;

	// The four discrete Laplacian sample points.
	double K_m = K - delta_K;
	double K_p = K + delta_K;
	double omega_m = omega - delta_omega;
	double omega_p = omega + delta_omega;

	// The time-summed laplacian
	double lap = 0.0;
	int nit = 0;

	for (int j=0; j<itermax/LAP_ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* First, we give a spin for 500 cycles, giving the non-chaotic
		 * parts a chance to phase-lock */
		for (int iter=0; iter<LAP_SETTLE_TIME; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
		}

		/* OK, now, we track five different points */
		double xom = x;
		double xop = x;
		double xkm = x;
		double xkp = x;
		for (int iter=0; iter < LAP_ITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
			xom += omega_m - K * sin (2.0 * M_PI * xom);
			xop += omega_p - K * sin (2.0 * M_PI * xop);
			xkm += omega - K_m * sin (2.0 * M_PI * xkm);
			xkp += omega - K_p * sin (2.0 * M_PI * xkp);

			// The sampling pattern
			lap += 4.0*x - xom - xop - xkm - xkp;
			nit ++;
		}
	}

	double lapavg = lap / ((double) nit);
	return lapavg;
}

/*-------------------------------------------------------------------*/
/*
 * Compute the metric distance between nearby orbits for the circle
 * map. This computes the l_1 or l_2 distance between the time series
 * for Lapalcians of points in the circle map.
 *
 * That is, consider an orbit (a time series) at fixed (K, omega). Then
 * consider the orbit for a nearby neighbor (K+deltaK, omega+deltaomega)
 * In mode-locked regions, we expect these to move together. At the
 * boundaries, we expect these to be uncorellated.
 *
 * At each point in the time-series, take the absolute value of the
 * distance between the series:
 *   dist = (1/N) sum_n abs(D_n(K, o))
 * where
 *   D_n(K, o) = x_n(K, o) - x_n(K+dK, o)
 *             + x_n(K, o) - x_n(K-dK, o)
 *             + x_n(K, o) - x_n(K, o+do)
 *             + x_n(K, o) - x_n(K, o-do)
 * The D_n above is the five-point Laplacian; more generally, the
 * (M+1) point Laplacian is used, with M=MET_SPOKES below.
 *
 * The above dist is the Banach l_1 distance aka Manhattan distance.
 * Alternately, the Euclidean (Hilbert) distance can be used:
 *   dist = sqrt( (1/N) sum_n |D_n(K, o)|^2 )
 * by compile-time adjustment; see below.
 */

#define MET_SETTLE_TIME 	0

// Iteration depth
// #define MET_ITER_DEPTH 120
// #define MET_ITER_DEPTH 480
#define MET_ITER_DEPTH 1920

// Neighborhood samples
#define MET_SPOKES 11

double
circle_metric (double omega, double K, int itermax, double param)

{
	// Sample offsets. We want this to be the distance to the neighboring
	// pixel, more or less. User-specified. A hard-coded 0.001 is not a
	// bad place to start.
	double delta = param; // 0.001;

	// The sampled regions.
	double Koff[MET_SPOKES];
	double omoff[MET_SPOKES];
	double t = rand();
	t /= RAND_MAX;
	for (int k=0; k<MET_SPOKES; k++)
	{
		Koff[k] = K + delta*cos(2.0 * M_PI * (k+t) / ((double) MET_SPOKES));
		omoff[k] = omega + delta*sin(2.0 * M_PI * (k+t) / ((double) MET_SPOKES));
	}

	// The time-summed distance
	double dist = 0.0;
	int nit = 0;

	for (int j=0; j<itermax/MET_ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* First, we give a spin for 500 cycles, giving the non-chaotic
		 * parts a chance to phase-lock */
		for (int iter=0; iter<MET_SETTLE_TIME; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
		}

		/* OK, now, we track MET_SPOKES+1 different points */
		double xoff[MET_SPOKES];
		for (int k=0; k<MET_SPOKES; k++)
		{
			xoff[k] = x;
		}

		for (int iter=0; iter < MET_ITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);

			// Lappy is the averaged Laplacian
			double lappy = 0.0;
			for (int k=0; k<MET_SPOKES; k++)
			{
				xoff[k] += omoff[k] - Koff[k] * sin (2.0 * M_PI * xoff[k]);
				lappy += x-xoff[k];
			}
			lappy /= (double) MET_SPOKES;

// #define MET_L1_METRIC 1
#ifdef MET_L1_METRIC
			dist += fabs(lappy);
#endif
#define MET_L2_METRIC 1
#ifdef MET_L2_METRIC
			dist += lappy*lappy;
#endif
			nit ++;
		}
	}

	double met = dist / ((double) nit);
#ifdef MET_L2_METRIC
	met = sqrt(met);
#endif
	return met;
}

/*-------------------------------------------------------------------*/
/*
 * Compute an averaged gradient between nearby orbits for the circle
 * map. This computes the l_1 or l_2 distance between the time series
 * for two nearby points in the circle map.
 *
 * As above, but instead an averaged gradient is used.
 *
 * At each point in the time-series, take the gradient between two
 * neighboring points in the series:
 *   G_n(K, o) = x_n(K, o) - x_n(K+dK, o+do)
 * To avoid sampling bias with respect to the direction, an average is
 * taken:
 *   GA_n(K, o) = (1/M) sum_{dK,do} |G_n(K, o)|
 * for M different directions dK,do.  Note the abssolute value above:
 * this becomes an averaged gradient, behaving like a first derivative;
 * without the absolute value, this would behave like a second
 * derivative, i.e. a Laplacian.
 *
 * The rest proceeds as before, so that
 *   dist = (1/N) sum_n abs(GA_n(K, o))
 * is the Banach l_1 distance, etc.
 */

#define GRD_SETTLE_TIME 	0

// Iteration depth
// #define GRD_ITER_DEPTH 120
// #define GRD_ITER_DEPTH 480
#define GRD_ITER_DEPTH 1920

// Neighborhood samples
#define GRD_SPOKES 11

double
circle_gradient (double omega, double K, int itermax, double param)

{
	// Sample offsets. We want this to be the distance to the neighboring
	// pixel, more or less. User-specified. A hard-coded 0.001 is not a
	// bad place to start.
	double delta = param; // 0.001;

	// The sampled regions.
	double Koff[GRD_SPOKES];
	double omoff[GRD_SPOKES];
	double t = rand();
	t /= RAND_MAX;
	for (int k=0; k<GRD_SPOKES; k++)
	{
		Koff[k] = K + delta*cos(2.0 * M_PI * (k+t) / ((double) GRD_SPOKES));
		omoff[k] = omega + delta*sin(2.0 * M_PI * (k+t) / ((double) GRD_SPOKES));
	}

	// The time-summed distance
	double dist = 0.0;
	int nit = 0;

	for (int j=0; j<itermax/GRD_ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* First, we give a spin for 500 cycles, giving the non-chaotic
		 * parts a chance to phase-lock */
		for (int iter=0; iter<GRD_SETTLE_TIME; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
		}

		/* OK, now, we track GRD_SPOKES+1 different points */
		double xoff[GRD_SPOKES];
		for (int k=0; k<GRD_SPOKES; k++)
		{
			xoff[k] = x;
		}

		for (int iter=0; iter < GRD_ITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
			double grad = 0.0;
			for (int k=0; k<GRD_SPOKES; k++)
			{
				xoff[k] += omoff[k] - Koff[k] * sin (2.0 * M_PI * xoff[k]);
#define GRD_L1_METRIC 1
#ifdef GRD_L1_METRIC
				grad += fabs(x-xoff[k]);
#endif
// #define GRD_L2_METRIC 1
#ifdef GRD_L2_METRIC
				grad += (x-xoff[k]) * (x-xoff[k]);
#endif
			}
			grad /= (double) GRD_SPOKES;
#ifdef GRD_L2_METRIC
			grad = sqrt(grad);
#endif
			nit ++;
#ifdef VEC_L1_METRIC
			dist += grad;
#endif
#define VEC_L2_METRIC
#ifdef VEC_L2_METRIC
			dist += grad*grad;
#endif
		}
	}

	double met = dist / ((double) nit);
#ifdef VEC_L2_METRIC
	met = sqrt(met);
#endif
	return met;
}

/*-------------------------------------------------------------------*/
/*
 * Much like above, but we compare orbits for +K to -K.
 */

#define FLP_SETTLE_TIME 	0

// Iteration depth
#define FLP_ITER_DEPTH 120
// #define FLP_ITER_DEPTH 480
// #define FLP_ITER_DEPTH 1920

double
circle_flip (double omega, double K, int itermax, double param)

{
	// The time-summed distance
	double dist = 0.0;
	int nit = 0;

	for (int j=0; j<itermax/FLP_ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;
		double y = t;

		/* First, we give a spin for 500 cycles, giving the non-chaotic
		 * parts a chance to phase-lock */
		for (int iter=0; iter<FLP_SETTLE_TIME; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
			y += omega + K * sin (2.0 * M_PI * y);
		}

		for (int iter=0; iter < FLP_ITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
			y += omega + K * sin (2.0 * M_PI * y);
// #define FLP_L1_METRIC 1
#ifdef FLP_L1_METRIC
			dist += fabs(x-y);
#endif
#define FLP_L2_METRIC 1
#ifdef FLP_L2_METRIC
			dist += (x-y) * (x-y);
#endif
			nit ++;
		}
	}

	double met = dist / ((double) nit);
#ifdef FLP_L2_METRIC
	met = sqrt(met);
#endif
	return met;
}

/*-------------------------------------------------------------------*/
/* Bifurcation diagram callback, does one row at a time */

static void
bifurcation_diagram
(float *array,
	int array_size,
	double x_center,
	double x_width,
	double K,
	int itermax,
	double omega)
{
	double	x=0.0;
	int		iter,j;
	int cnt=0;

	/* clear out the row */
	for (j=0; j<array_size; j++)
	{
		array[j] = 0.0;
	}

#define BITER_DEPTH 500
	for (j=0; j<itermax/BITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		t -= 0.5;
		x = t;

		/* OK, now start iterating the circle map */
		for (iter=0; iter < BITER_DEPTH; iter++) {
			x += omega - K * sin (2.0 * M_PI * x);

			double en = array_size * (x-floor(x));
			int n = en;
			if (0 > n) n = 0;
			if (n >= array_size) n = array_size-1;
			array[n] += 1.0;
			cnt ++;
		}
	}

	for (j=0; j<array_size; j++)
	{
		array[j] /= cnt;
	}
}

/*-------------------------------------------------------------------*/

static double circle_map (double omega, double K, int itermax, double param)
{
	// return winding_number (omega, K, itermax);
	// return noisy_winding_number (omega, K, itermax, param);
	// return rms_winding_number (omega, K, itermax);
	// return circle_poincare_recurrance_time (omega, K, itermax);
	// return circle_laplacian (omega, K, itermax, param);
	return circle_gradient (omega, K, itermax, param);
	// return circle_metric (omega, K, itermax, param);
	// return circle_flip (omega, K, itermax, param);
}

DECL_MAKE_HEIGHT (circle_map);

// DECL_MAKE_BIFUR(bifurcation_diagram)

/* --------------------------- END OF LIFE ------------------------- */

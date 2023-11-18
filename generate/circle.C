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
	double	x=0.0;
	int		iter,j;
	int cnt=0;
	double start=0.0, end=0.0;

#define SAMP 150
	for (j=0; j<itermax/SAMP; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		x = t;
		start += x;

		/* OK, now start iterating the circle map */
		for (iter=0; iter < SAMP; iter++) {
			x += omega - K * sin (2.0 * M_PI * x);
			cnt ++;
		}
		end += x;
	}

	x = (end-start) / ((double) cnt);
	return x;
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

	for (j=0; j<itermax/SAMP; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		x = t;
		start += x;

		/* OK, now start iterating the circle map */
		for (iter=0; iter < SAMP; iter++) {
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
#define RSAMP 200

double
circle_poincare_recurrance_time (double omega, double K, int itermax)

{
	double	x, y;
	double	xpoint;
	int		j, iter;
	long		num_recurs, time_recur=0;

  	num_recurs = 0;
	for (j=0; j<itermax/RSAMP; j++)
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
		for (iter=0; iter < RSAMP; iter++)
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
 * Compute the laplacian of the orbit (the charge) for the circle map
 * This computes the series average for the four-point discrete
 * laplacian for neighboring orbits.
 *
 * That is, consider an orbit (a time series) at fixed (K, omega). Then
 * consider the orbit for a nearby neighbor (K+deltaK, omega+deltaomega)
 * In mode-locked regions, we expect these to move together. At the
 * boundaries, we expect these to be uncorellated.
 *
 * At each point in the time-series, take the four-point Laplacian, viz:
 * 4x_n(K, o) - x_n(K-dK, o) - x_n(K+dK, o) - x_n(K, o-do) - x_n(K, o+do)
 * and then just average this together over all n.
 */

#define LAP_SETTLE_TIME 	90
#define LAP_RSAMP 20

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

	for (int j=0; j<itermax/LAP_RSAMP; j++)
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
		for (int iter=0; iter < LAP_RSAMP; iter++)
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

#define BSAMP 500
	for (j=0; j<itermax/BSAMP; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		t -= 0.5;
		x = t;

		/* OK, now start iterating the circle map */
		for (iter=0; iter < BSAMP; iter++) {
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
	return circle_laplacian (omega, K, itermax, param);
}

DECL_MAKE_HEIGHT (circle_map);

// DECL_MAKE_BIFUR(bifurcation_diagram)

/* --------------------------- END OF LIFE ------------------------- */

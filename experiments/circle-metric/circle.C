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
 * more stuff -- November 2023
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
 * Simple, simplistic algo.
 * -- Iterate freely for SETTLE_TIME.
 * -- Record this as the starting point xpoint.
 * -- Count how long it takes to come back to within EPSILON of this
 *    starting point.
 * -- Average over multiple random starts.
 *
 * Problems with this algo: A long orbit might accidentally come
 * within epsilon of itself, and thus will be measured shorter than
 * it actually is. Seems like this might be especially a problem near
 * bifurcations and maybe at the transition to chaos?
 * The bin-counting algo below is meant to cure this issue.
 */

#define EPSILON  	0.002
// #define SETTLE_TIME 190
#define SETTLE_TIME 390
// #define SETTLE_TIME 1490
// #define RITER_DEPTH 500       // Iteration depth
#define RITER_DEPTH 3500
// #define RITER_DEPTH 18500

double
circle_poincare_recurrance_time (double omega, double K, int itermax)

{
	long num_recurs = 0;
	long time_recur=0;

	for (int j=0; j<itermax/RITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* First, we give a spin for 500 cycles, giving the non-chaotic
		 * parts a chance to phase-lock */
		for (int iter=0; iter<SETTLE_TIME; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
		}

		/* OK, now, we begin to measure the average amount of time
		 * to recur. Note that we don't have todo += with iter,
		 * since its already a running sum. */
		double xpoint = x;
		long ptime = 0;
		for (int iter=0; iter < RITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
			double y = fabs (x-xpoint);
			y -= floor (y);
			if (y < EPSILON)
			{
				num_recurs ++;
				ptime = iter;
			}
		}
		time_recur += ptime;
	}

	/* trec is the (normalized) number of cycles to reach recurrance */
	double trec = ((double) time_recur) / ((double) num_recurs);
	return trec;
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the circle map
 * Fancy bin-counting version.
 *
 * Algo:
 * -- iterate, placing into bins.
 * -- count number of non-empty bins
 * -- recount with min threshold
 *
 * Seems like it should be able to avoid some of the issues with
 * the naive recurrence-time algo, above. Of course, it does have
 * its own issues, e.g. with noise, jitter.
 */

// Make sure that ITER_DEPTH is at least 4x larger than nbins,
// so that for chaotic orbits, each bin is hit maybe 4 times.
// This will avoid issues with the thresholding, below.
// #define PNC_NBINS 400
// #define PNC_NBINS 800
#define PNC_NBINS 5000
#define PNC_SETTLE_TIME 190
// #define PNC_ITER_DEPTH 1920       // Iteration depth
// #define PNC_ITER_DEPTH 8111
#define PNC_ITER_DEPTH 58111

double
circle_poincare_bincount (double omega, double K, int itermax)

{
	long period_len = 0;
	long nmeasurements = 0;

	for (int j=0; j<itermax/PNC_ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* First, we give a spin for 500 cycles, giving the non-chaotic
		 * parts a chance to phase-lock */
		for (int iter=0; iter<PNC_SETTLE_TIME; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
		}

		double bins[PNC_NBINS+1];
		for (int ib=0; ib<=PNC_NBINS; ib++)
			bins[ib] = 0;

		// OK, now, bin-count as we move along.
		// bin-counting is always modulo one, because that is
		// all that sine cares about.
		for (int iter=0; iter < PNC_ITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
			double yb = PNC_NBINS * (x - floor (x));
			int ib = (int) yb;
			bins[ib]++;
		}

		// Count number of non-empty bins.
		int nhits = 0;
		for (int ib=0; ib<=PNC_NBINS; ib++)
			if (0 < bins[ib]) nhits++;

		// Recount, using a threshold.
		double binavg = ((double) PNC_ITER_DEPTH) / ((double) nhits);
		int ibavg = (int) (0.5*binavg);

		// Final count, rejecting mostly empty bins.
		int cyclen = 0;
		for (int ib=0; ib<=PNC_NBINS; ib++)
			if (ibavg < bins[ib]) cyclen++;

		// Tally up.
		period_len += cyclen;
		nmeasurements ++;
	}

	/* plen is the average length of the period cycle. */
	double plen = ((double) period_len) / ((double) nmeasurements);
	return plen;
}

/*-------------------------------------------------------------------*/
/*
 * Compute the support of the measure for the circle map. For
 * periodic orbits, this provides a number that is more or less the
 * same as the poincare recurrence time (divided by the number of bins).
 * For the chaotic orbits, it returns an approximation to the support
 * of the invariant measure, as a fraction running from zero to one.
 *
 * Algo:
 * -- iterate, placing into bins.
 * -- count number of non-empty bins
 */

// Make sure that ITER_DEPTH is at least 4x larger than nbins,
// so that for chaotic orbits, each bin is hit maybe 4 times.
// This will avoid issues with the thresholding, below.
#define ISS_NBINS 1000
#define ISS_SETTLE_TIME 190
// #define ISS_ITER_DEPTH 1920       // Iteration depth
#define ISS_ITER_DEPTH 8111
// #define ISS_ITER_DEPTH 38111
// #define ISS_THRESH 0.05  // reject threshold for empty bins
#define ISS_THRESH 0.03

double
circle_support (double omega, double K, int itermax)
{
	long non_empty_bins = 0;
	long nmeasurements = 0;

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

		double bins[ISS_NBINS+1];
		for (int ib=0; ib<=ISS_NBINS; ib++)
			bins[ib] = 0;

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

		// Count number of non-empty bins. Reject almost-empty bins.
		double binavg = ((double) ISS_ITER_DEPTH) / ((double) ISS_NBINS);
		int ithresh = (int) floor (ISS_THRESH*binavg);

		// Count rejecting mostly empty bins.
		int support = 0;
		for (int ib=0; ib<=ISS_NBINS; ib++)
			if (ithresh < bins[ib]) support++;

		// Tally up.
		non_empty_bins += support;
		nmeasurements ++;
	}

	// supp is the average number of non-empty bins
	double supp = ((double) non_empty_bins) / ((double) nmeasurements);

	// Divide by the total number of bins, to get the support of the
	// measure. The support goes from zero to one.
	supp /= (double) ISS_NBINS;
	return supp;
}

/*-------------------------------------------------------------------*/
/*
 * Perform an estimate of the recurrence time as the reciprocal of
 * the measure in a given interval.
 *
 * Algo:
 * -- iterate, approximating the measure via bin-count.
 * -- return content of a specific bin.
 */

// Make sure that ITER_DEPTH is at least 4x larger than nbins,
// so that for chaotic orbits, each bin is hit maybe 4 times.
// This will avoid issues with the thresholding, below.
#define MEA_NBINS 500
#define MEA_SETTLE_TIME 190
// #define MEA_ITER_DEPTH 1920       // Iteration depth
#define MEA_ITER_DEPTH 8111
// #define MEA_ITER_DEPTH 38111

double
circle_poincare_measure(double omega, double K, int itermax)
{
	double bins[MEA_NBINS+1];
	for (int ib=0; ib<=MEA_NBINS; ib++)
		bins[ib] = 0;

	for (int j=0; j<itermax/MEA_ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* First, we give a spin for 500 cycles, giving the non-chaotic
		 * parts a chance to phase-lock */
		for (int iter=0; iter<MEA_SETTLE_TIME; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
		}

		// OK, now, bin-count as we move along.
		// bin-counting is always modulo one, because that is
		// all that sine cares about.
		for (int iter=0; iter < MEA_ITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
			double yb = MEA_NBINS * (x - floor (x));
			int ib = (int) yb;
			bins[ib]++;
		}
	}

	// Count grand total
	long totcnt = 0;
	for (int ib=0; ib<=MEA_NBINS; ib++)
		totcnt += bins[ib];

	double totmeas = 0.0;
	int nsamp = 0;

	// Pick bins at random, with probability given by the measure.
	for (int j=0; j<itermax/MEA_ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		// Thermalize
		for (int iter=0; iter<MEA_SETTLE_TIME; iter++)
			x += omega - K * sin (2.0 * M_PI * x);

		// Report
		for (int iter=0; iter < MEA_ITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
			double yb = MEA_NBINS * (x - floor (x));
			int ib = (int) yb;

			if (0 < bins[ib])
			{
				double meas = bins[ib];
				meas /= (double) totcnt;
				meas *= (double) MEA_NBINS;
				totmeas += meas;
				nsamp ++;
			}
		}
	}

	double avgmeas = totmeas / ((double) nsamp);

	// Recurrence time estimate
	double rtime = 1.0 / avgmeas;
	return rtime;
}

/*-------------------------------------------------------------------*/
/*
 * Compute the Lyapunov exponent for the circle map. Cheap and dirty.
 */

#define LYA_SETTLE_TIME 90  // Settle time
// #define LYA_SETTLE_TIME 490  // settle time
#define LYA_ITER_DEPTH 1920       // Iteration depth
// #define LYA_ITER_DEPTH 19311

double
circle_lyapunov (double omega, double K, int itermax)
{
	double totlya = 0.0;
	int nobs = 0;

	for (int j=0; j<itermax/LYA_ITER_DEPTH; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		double x = t;

		/* First, we give a spin for 500 cycles, giving the non-chaotic
		 * parts a chance to phase-lock */
		for (int iter=0; iter<LYA_SETTLE_TIME; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);
		}

		for (int iter=0; iter < LYA_ITER_DEPTH; iter++)
		{
			x += omega - K * sin (2.0 * M_PI * x);

			// "push" is the differential push contributing
			// to the non-linear changes in the orbit. Take
			// before iterating; use the bin after iterating.
			// (So that it's correctly convolved w/measure.)
			double push = abs(1.0 - 2.0 * M_PI * K * cos (2.0 * M_PI * x));
			totlya += push;
			nobs ++;
		}
	}

	double avglya = totlya / nobs;
	avglya = log(avglya);

	// That's intersting. Lets divide by K
	if (0.0 != K) avglya /= K;
	// if (0.0 != K) avglya /= pow(K, 1.5);
	// if (0.0 != K) avglya /= log(1+K);
	// if (0.0 != K) avglya -= 0.25*log(K);

	return avglya;
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
#define VEC_L1_METRIC
#ifdef VEC_L1_METRIC
			dist += grad;
#endif
// #define VEC_L2_METRIC
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

static double circle_map (double a, double b, int itermax, double param)
{
#define USUAL 1
#ifdef USUAL
	double omega = a;
	double K = b;
#else
	// Spin by 45 degrees.
	double omega = a-b;
	double K = a+b;
#endif

	// return winding_number (omega, K, itermax);
	// return noisy_winding_number (omega, K, itermax, param);
	// return rms_winding_number (omega, K, itermax);
	// return circle_poincare_recurrance_time (omega, K, itermax);
	// return circle_poincare_bincount (omega, K, itermax);
	// return circle_poincare_measure (omega, K, itermax);
	// return circle_support (omega, K, itermax);
	return circle_lyapunov (omega, K, itermax);
	// return circle_laplacian (omega, K, itermax, param);
	// return circle_gradient (omega, K, itermax, param);
	// return circle_metric (omega, K, itermax, param);
	// return circle_flip (omega, K, itermax, param);
}

DECL_MAKE_HEIGHT (circle_map);

// DECL_MAKE_BIFUR(bifurcation_diagram)

/* --------------------------- END OF LIFE ------------------------- */

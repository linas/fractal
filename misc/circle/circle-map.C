/*
 * circle-map.C
 *
 * Brand-new recreation of the circle-map stuff.
 * The original version is in generate/circle.C
 * Dec 2017
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

// Traditional circle map. Kbar = K/2pi
double circle_map(double xn, double omega, double Kbar)
{
	double xnp1 = xn + omega - Kbar * sin(2.0 * M_PI * xn);
	return xnp1;
}

// Triangle-map approximation to circle map. Kbar = K/2pi
double triangle_map(double xn, double omega, double Kbar)
{
	double tri = xn;
	tri -= floor(tri);

	if (0.75 < tri) tri = tri - 1.0;
	else if (0.25 < tri) tri = 0.5 - tri;

	double K = Kbar * 4.0;
	double xnp1 = xn + omega - K * tri;
	return xnp1;
}

/*-------------------------------------------------------------------*/
/*
 * This routine computes average winding number taken by the map.
 */
double winding_number(double omega, double Kbar, int itermax,
                      double (func)(double, double, double) )
{
   double	x=0.0;
   int		iter,j;
	int cnt=0;
	double start=0.0, end=0.0;

#define SAMP 750
	for (j=0; j<itermax; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		x = t;
		start += x;

		/* OK, now start iterating the circle map */
		for (iter=0; iter < SAMP; iter++) {
			x = func(x, omega, Kbar);
			cnt ++;
		}
		end += x;
	}

	x = (end-start) / ((double) cnt);
	return x;
}

/*-------------------------------------------------------------------*/
/*
 * Compute the poincare recurrance time for the circle map
 */

// #define EPSILON 0.001
#define EPSILON 0.003
#define SETTLE_TIME 190
#define RSAMP 2400

double
recurrance_time (double omega, double Kbar, int itermax,
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


static double circle_gram(double omega, double Kbar, int itermax, double param)
{
	// return winding_number(omega, Kbar, itermax, circle_map);
	// return winding_number(omega, Kbar, itermax, triangle_map);
	// return recurrance_time(omega, Kbar, itermax, circle_map);
	return recurrance_time(omega, Kbar, itermax, triangle_map);
}

DECL_MAKE_HEIGHT (circle_gram);

#if 0
int main ( int argc, char * argv[])
{
	double om = atof(argv[1]);
	double kb = atof(argv[2]);
	winding_number(om, kb, 150, circle_map);
}
#endif

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
double sawtooth_map(double xn, double omega, double Kbar)
{
	double tri = xn;
	while (tri < 0.0) tri += 1.0;
	while (tri > 1.0) tri -= 1.0;

	if (0.75 < xn) tri = xn - 1.0;
	else if (0.25 < xn) tri = 0.5 - xn;

	double xnp1 = xn + omega + Kbar * tri;
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

#define SAMP 150
	for (j=0; j<itermax/SAMP; j++)
	{
		double t = rand();
		t /= RAND_MAX;
		x = t;
		start += x - 0.5; // subtract 1/2 so that avg is zero.

		/* OK, now start iterating the circle map */
		for (iter=0; iter < SAMP; iter++) {
			x = func(omega, Kbar, x);
			cnt ++;
		}
		end += x;
	}

	x = (end-start) / ((double) cnt);
printf("duude windw= %g for  %g %g \n", x, omega, Kbar);
	return x;
}


static double circle_gram(double omega, double Kbar, int itermax, double param)
{
printf("duuude %g %g \n", omega, Kbar);
	return winding_number(omega, Kbar, itermax, circle_map);
}

DECL_MAKE_HEIGHT (circle_gram);

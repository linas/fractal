/*
 * circle-map.C
 *
 * Brand-new recreation of the circle-map stuff.
 * Dec 2017
 */

#include <math.h>

// Traditional circle map. Kbar = K/2pi
double circle_map(double xn, double omega, double Kbar)
{
	double xnp1 = xn + omega - Kbar * sin(2.0 * M_PI * xn);
	while (xnp1 < 0.0) xnp1 += 1.0;
	while (xnp1 > 1.0) xnp1 -= 1.0;
	return xnp1;
}

// Triangle-map approximation to circle map. Kbar = K/2pi
double sawtooth_map(double xn, double omega, double Kbar)
{
	double tri = xn;
	if (0.75 < xn) tri = xn - 1.0;
	else if (0.25 < xn) tri = 0.5 - xn;

	double xnp1 = xn + omega + Kbar * tri;
	while (xnp1 < 0.0) xnp1 += 1.0;
	while (xnp1 > 1.0) xnp1 -= 1.0;
	return xnp1;
}

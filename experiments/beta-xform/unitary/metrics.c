/*
 * metrics.c
 * Distances between pairs of points in various hyperbolic geometry
 * models. Poincare Disk, Klien geometry disk.
 *
 * February 2024
 */

#include <math.h>

// Isometric invariant on the Poincare disk
// See https://en.wikipedia.org/wiki/Poincar%C3%A9_disk_model
// It is delta(u,v) = 2 |u-v|^2 / (1-|u|^2)(1-|v|^2)
// Input args are modulus and phase normalized to unit interval.
double poincare_iso(double ru, double tu, double rv, double tv)
{
	double ux = ru * cos(tu * 2.0*M_PI);
	double uy = ru * sin(tu * 2.0*M_PI);
	double vx = rv * cos(tv * 2.0*M_PI);
	double vy = rv * sin(tv * 2.0*M_PI);

	double delta = 2.0 * ((ux-vx)*(ux-vx) + (uy-vy)*(uy-vy));
	delta /= 1.0 - (ux*ux + uy*uy);
	delta /= 1.0 - (vx*vx + vy*vy);

	return delta;
}

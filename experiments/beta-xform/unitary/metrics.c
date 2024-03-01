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

	double up = 1.0 - (ux*ux + uy*uy);
	double vp = 1.0 - (vx*vx + vy*vy);

	if (up < 1.0e-10)
		fprintf(stderr, "point u on the boundary: r=%f phi=%f\n", ru, tu);
	if (vp < 1.0e-10)
		fprintf(stderr, "point v on the boundary: r=%f phi=%f\n", rv, tv);

	double delta = 2.0 * ((ux-vx)*(ux-vx) + (uy-vy)*(uy-vy));
	delta /= up * vp;

	return delta;
}

// Assume inputs are in Klein coords
double klein_iso(double kru, double tu, double krv, double tv)
{
	double ru = kru / (1.0 + sqrt(1.0 - kru*kru));
	double rv = krv / (1.0 + sqrt(1.0 - krv*krv));
	return poincare_iso(ru, tu, rv, tv);
}

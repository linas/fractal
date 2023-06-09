/*
 * coord-xforms.C
 *
 * FUNCTION:
 * provide coordinate transforms from upper half plane to the 
 * poincare disk and the q-series disk, and thier inverses.
 *
 * HISTORY:
 * New, May 2005
 */

#include <math.h>
#include "coord-xforms.h"

/**
 * poincare_disk_to_plane_coords
 *
 * Return i (1+z)/(1-z) where z=x+iy is a disk coord
 * Interior of disk is mapped only to the upper half-plane,
 * and the real number line is wrapped around the circumference.
 * Exterior of disk maps to the lower half-plane.
 */
void 
poincare_disk_to_plane_coords (double x, double y, 
                               double *px, double *py)
{
	double deno = (1.0-x)*(1.0-x) + y*y;
	deno = 1.0/deno;
	*px = -2.0*y*deno;
	*py = (1.0 - x*x - y*y) * deno;
}

/**
 * plane_to_poincare_disk_coords
 *
 * Return (tau-i)/(tau+i) where tau=x+iy are (upper half) plane coords
 * Points in the upper half plane map to inside the unit circle, 
 * Points in teh lower-half plane map to the outside of the circle.
 */
void 
plane_to_poincare_disk_coords (double x, double y, 
                               double *px, double *py)
{
	double deno = x*x + (y+1.0)*(y+1.0);
	deno = 1.0/deno;
	*px = (x*x + y*y - 1.0) * deno;
	*py = -2.0*x*deno; 
}

void 
plane_to_q_disk_coords (double tau_re, double tau_im, 
                        double *px, double *py)
{
	/* now go back to q-series coords */
	double rq = exp (-tau_im * 2.0 * M_PI);
	*px = rq * cos (tau_re * 2.0 * M_PI);
	*py = rq * sin (tau_re * 2.0 * M_PI);
}

void 
q_disk_to_plane_coords (double qre, double qim,
                        double *px, double *py)
{
	/* map from q-series coords to the upper half-plane */
	double tau_im = -log (sqrt (qre*qre +qim*qim)) / (2.0*M_PI);
	double tau_re = atan2 (qim, qre) /(2.0*M_PI);
	*px = tau_re;
	*py = tau_im;
}

void 
mobius_xform (double a, double b, double c, double d,
              double tau_re, double tau_im, 
              double *px, double *py)
{
	/* now apply mobius */
	double deno = c*tau_re+d;
	deno = deno*deno + c*c*tau_im*tau_im;
	deno = 1.0/deno;

	tau_re = (a*tau_re+b)*(c*tau_re+d) + a*c*tau_im*tau_im;
	*px = tau_re * deno;
	*py = (a*d-b*c) * tau_im * deno;
}


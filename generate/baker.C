
/*
 * baker.C
 *
 * Show invariants and eigenvectors of Baker's map.
 *
 * Linas Vepstas Sept 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

/* Ising model nerest-neighbor interaction */
double ising(double x, double y)
{
	/* higher energy if pair of spins is  aligned, else, energy is lower */
	if ((x<0.25) || (x>=0.75)) return 0.3;
	return -0.3;
}

double parabola(double x, double y)
{
	x -= 0.5;
	y -= 0.5;
	return 1.5*(x*x+y*y-0.125);
}

double b1(double x, double y)
{
	return y-0.5;
}

/* =============================================================== */

static int niter=10; 
static double tee = 0.3;

void eigenvec (double *re, double *im, double (*func)(double, double), 
                double x, double y)
{
	int i;
	double hre = 0.0;
	double him = 0.0;
	
	double xo = x;
	double yo = y;
	double w = 1.0;
#if 0
	for (i=0; i<niter+1; i++)
	{
		// w = 1.0/(i-tee);
		hre += w*func (x,y);
		int ty = (int) floor(2.0*y);
		x = 0.5*(x+ty);
		y = 2.0*y-ty;

	}
#endif

	x = yo;
	y = xo;
	for (i=0; i<niter; i++)
	{
		int ty = (int) floor(2.0*y);
		x = 0.5*(x+ty);
		y = 2.0*y-ty;
		// w = 1.0/(-(i+1)-tee);
		hre += w*func (y,x);
		w *= tee;
	}

	*re = hre;
	*im = him;
}

static double 
density (double x, double y, int itermax, double param)
{
	double re, im;

	niter = itermax;
	tee = param;
	eigenvec (&re, &im, ising, x, y);
	// eigenvec (&re, &im, parabola, x, y);
	// eigenvec (&re, &im, b1, x, y);
	re = exp (-re);
	return re;
}

DECL_MAKE_HEIGHT(density);

/* --------------------------- END OF LIFE ------------------------- */

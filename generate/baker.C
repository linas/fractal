
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

double ising(double x, double y)
{
	// x = y;
	if ((x<0.25) || (x>=0.75)) return 0.3;
	return -0.3;
}

void eigenvec (double *re, double *im, double (*func)(double, double), 
                double x, double y)
{
	int i;
	double hre = 0.0;
	double him = 0.0;
	
	double xo = x;
	double yo = y;
#define NITER 10
	for (i=0; i<NITER; i++)
	{
		hre += func (x,y);
		int ty = (int) floor(2.0*y);
		x = 0.5*(x+ty);
		y = 2.0*y-ty;
	}

	x = yo;
	y = xo;
	for (i=0; i<NITER; i++)
	{
		int ty = (int) floor(2.0*y);
		x = 0.5*(x+ty);
		y = 2.0*y-ty;
		hre += func (x,y);
	}
	// hre /= 80;
	// him /= 80;

	*re = hre;
	*im = him;
}

static double 
density (double x, double y, int itermax, double param)
{
	double re, im;

	eigenvec (&re, &im, ising, x, y);
	re = exp (-re);
	return re;
}

DECL_MAKE_HEIGHT(density);

/* --------------------------- END OF LIFE ------------------------- */

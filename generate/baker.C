
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
	if ((x<0.25 || (x>0.75)) return -0.3;
	return 0.3;
}

void eigenvec (double *re, double *im, double (*func)(double, double), 
                double x, double y)
{
	int i;
	double hre = him = 0.0;
	
	for (i=0; i<40; i++)
	{
		hre += func (x,y);
		int ty = 2.0*y;
		x = 0.5*(x+ty);
		y = 2.0*y-ty;
	}

	*re = hre;
}

static double 
density (double x, double y, int itermax, double param)
{
	double re, im;

	eigenvec (&re, &im, ising, x, y);
	return re;
}

DECL_MAKE_HEIGHT(density);

/* --------------------------- END OF LIFE ------------------------- */

/*
 * gauss-red.C
 *
 * Perform Gaussian reduction i.e. iterate on
 * U(z) = 1/z - int Re(1/z) 
 *
 * Linas Vepstas January 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

static double gauss_red (double x, double y, int itermax, double param)
{
	int cnt = 0;

	double rad = (x-0.5)*(x-0.5) + y*y;
	while (rad < 0.25)
	{
		double mod = 1.0 / (x*x + y*y);
		x *= mod;
		y = -y * mod;

		x -= floor(x);

		rad = (x-0.5)*(x-0.5) + y*y;
		cnt ++;

		if (cnt > itermax) break;
	}

	return cnt;
}

DECL_MAKE_HEIGHT(gauss_red);

/* --------------------------- END OF LIFE ------------------------- */


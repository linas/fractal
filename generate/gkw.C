/*
 * gkw.C
 *
 * FUNCTION:
 * Display GKW operator, one matrix element per pixel.
 *
 * HISTORY:
 * New Jan 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../misc/continued-frobenius/ache.h"
#include "brat.h"


static double gkw_operator (double x, double y, int itermax, double param)
{
	int m = 100.0 * x;
	int p = 100.0 * y;
	double gkw = ache_mp(m,p);
	return gkw;
}

DECL_MAKE_HEIGHT(gkw_operator);

/* --------------------------- END OF LIFE ------------------------- */

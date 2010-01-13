/*
 * gkw.C
 *
 * FUNCTION:
 * Display GKW operator, one matrix element per pixel.
 *
 * HISTORY:
 * December 2003
 * Pixelize Jan 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "binomial.h"
#include "harmonic.h"

// Return the matrix element for H_mp aka the matrix element of GKW.
//
long double
ache_mp(int m, int p)
{
	int k;

	long double acc = 0.0L;
	long double sign = 1.0L;
	for (k=0; k<=p; k++)
	{
		long double term = zetam1 (k+m+2);
		term *= binomial (m+k+1,m);
		term *= binomial (p,k);
		term *= sign;
		acc += term;
		sign = -sign;
	}
	return acc;
}


static double gkw_operator (double x, double y, int itermax, double param)
{
	int p = 100.0 * x + 0.5;
	int m = 100.0 * y + 0.5;
	m = 100 - m;
	double gkw = ache_mp(m,p);
	gkw = fabs(gkw);
	gkw = log(gkw);
	gkw += 15.0;
	gkw /= 1.0;
// printf ("%d %d %f\n", m, p, gkw);
	return gkw;
}

DECL_MAKE_HEIGHT(gkw_operator);

/* --------------------------- END OF LIFE ------------------------- */

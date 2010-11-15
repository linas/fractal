/*
 * haar.C
 *
 * Har wavelets on complex plane
 *
 * Linas Vepstas November 2010
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "brat.h"

/* ======================================================================= */

	
static double haar (double re_q, double im_q, int itermax, double param)
{
	int n;

	double complex z = re_q + I*im_q;

	double complex acc = 0.0;
	for (n=0; n<4000; n++)
	{
		double complex term = csin(2.0*M_PI*(2*n+1)*z);
		term /= (double) (2*n+1);
		acc += term; 
	}

	double phase = carg(acc);

	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

DECL_MAKE_HEIGHT(haar);

/* --------------------------- END OF LIFE ------------------------- */

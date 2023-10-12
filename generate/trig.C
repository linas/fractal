/*
 * trig.C
 *
 * FUNCTION:
 * display trignometric functions
 *
 * HISTORY:
 * quick hack -- Linas Vepstas October 1989
 * modernize -- Linas Vepstas March 1996
 * more stuff -- January 2000
 * more stuff -- October 2004
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

static void trig_series_c (double re_q, double im_q, double *prep, double *pimp)
{
	double tmp;

	double rcz = cos (re_q) * cosh (im_q);
	double icz = sin (re_q) * sinh (im_q);

	*prep = rcz;
	*pimp = icz;
}

static double trig_series (double re_q, double im_q, int itermax, double param)
{
	double rep, imp;
	trig_series_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	// return rep;
	return (atan2 (imp,rep)+M_PI)/(2.0*M_PI);
}

DECL_MAKE_HEIGHT(trig_series);

/* --------------------------- END OF LIFE ------------------------- */

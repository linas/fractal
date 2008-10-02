/*
 * gkd.C
 *
 * FUNCTION:
 * Display generating function of the Gauss-Kuzmin distribution
 * (e.g. the moment generating function, etc.)
 *
 * HISTORY:
 * Linas - October 2008
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

static void gkd_series_c (double re_q, double im_q, double *prep, double *pimp)
{
	int k;
	
	double r = sqrt(re_q*re_q + im_q*im_q);
	double theta = atan2 (im_q,re_q);

	double rk = 1.0;
	double tk = 0.0;

	double re_g = 0.0;
	double im_g = 0.0;
	for (k=1; k<5000; k++)
	{
		double lg = 1.0 - 1.0 / ((k+1.0)*(k+1.0));
		lg = - log(lg);
		re_g = rk * cos(tk) * lg;
		re_g = rk * sin(tk) * lg;

		rk *= r;
		tk += theta;
	}

	*prep = re_g;
	*pimp = im_g;
}

static double gkd_series (double re_q, double im_q, int itermax, double param)
{
	double rep, imp;
	gkd_series_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	// return rep;
	return (atan2 (imp,rep)+M_PI)/(2.0*M_PI);
}

DECL_MAKE_HEIGHT(gkd_series);

/* --------------------------- END OF LIFE ------------------------- */

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

static void gkd_generic (double re_q, double im_q, double *prep, double *pimp)
{
	int k;

	double r = re_q;
	double theta = im_q;

	double rk = r;
	double cs = cos(theta);
	double sn = sin(theta);
	double ck = cs;
	double sk = sn;

	double re_g = 0.0;
	double im_g = 0.0;
	for (k=1; k<5000; k++)
	{
		double lg = 1.0 - 1.0 / ((k+1.0)*(k+1.0));
		lg = - log(lg);
		re_g += rk * ck * lg;
		im_g += rk * sk * lg;

		rk *= r;
		double tmp = ck*sn + sk*cs;
		ck = ck*cs - sk*sn;
		sk = tmp;
	}
	re_g /= M_LN2;
	im_g /= M_LN2;

	*prep = re_g;
	*pimp = im_g;
}

// Probability generating function
static void gkd_prob_gen_c (double re_q, double im_q, double *prep, double *pimp)
{
	double r = sqrt(re_q*re_q + im_q*im_q);
	double theta = atan2 (im_q,re_q);

	gkd_generic (r, theta, prep, pimp);
}

// Moment generating function
static void gkd_mom_gen_c (double re_q, double im_q, double *prep, double *pimp)
{
	double r = exp(re_q);
	double theta = im_q;

	gkd_generic (r, theta, prep, pimp);
}


static double gkd_series (double re_q, double im_q, int itermax, double param)
{
	double rep, imp;
	rep = re_q;
	imp = im_q;
	// gkd_mom_gen_c (re_q, im_q, &rep, &imp);
	// gkd_prob_gen_c (re_q, im_q, &rep, &imp);
	
	return sqrt (rep*rep+imp*imp);
	// return rep;
	// return (atan2 (imp,rep)+M_PI)/(2.0*M_PI);
}

DECL_MAKE_HEIGHT(gkd_series);

/* --------------------------- END OF LIFE ------------------------- */

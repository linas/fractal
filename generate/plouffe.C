/*
 * plouffe.C
 *
 * FUNCTION:
 * display plouffe-type q-series
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


static void totient_series_c (double re_q, double im_q, double *prep, double *pimp)
{
	double tmp;
	int i;
	*prep = 0.0;
	*pimp = 0.0;

	double rep = 0.0;
	double imp = 0.0;

	double qpr = 1.0;
	double qpi = 0.0;

	double qpmod = re_q*re_q+im_q*im_q;
	if (1.0 <= qpmod) return;

	for (i=0; i<max_terms; i++)
	{
		double t = totient_phi (i+1);

		t *= (i+1);
		t *= (i+1);
		t *= (i+1);

		rep += qpr *t;
		imp += qpi *t;

		/* compute q^k */
		tmp = qpr*re_q - qpi * im_q;
		qpi = qpr*im_q + qpi * re_q;
		qpr = tmp;

		qpmod = qpr*qpr + qpi*qpi;
		if (qpmod < 1.0e-30) break;
	}
	if (max_terms-1 < i)
	{
		// printf ("not converged re=%g im=%g modulus=%g\n", re_q, im_q, qpmod);
	}

#if NO_NOT_THIS
	/* multiply by (1-q)^2 */
	qpr = 1.0 - re_q;
	qpi = - im_q;
	tmp = qpr*qpr - qpi * qpi;
	qpi = 2.0*qpr*qpi;
	qpr = tmp;

	tmp = qpr*rep - qpi * imp;
	imp = qpr*imp + qpi * rep;
	rep = tmp;
#endif
	/* multiply by (1-|q|)^2 */
	tmp = 1.0 - sqrt (re_q*re_q + im_q*im_q);
	tmp *= tmp;

	rep *= tmp;
	imp *= tmp;

	*prep = rep;
	*pimp = imp;
}

static double totient_series (double re_q, double im_q, int itermax, double param)
{
	max_terms = itermax;
	double rep, imp;
	totient_series_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	return rep;
}

DECL_MAKE_HISTO(totient_series);

/* --------------------------- END OF LIFE ------------------------- */

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

static double ess;

static void plouffe_series_c (double re_q, double im_q, double *prep, double *pimp)
{
	double tmp;
	int n;
	*prep = 0.0;
	*pimp = 0.0;

	double rep = 0.0;
	double imp = 0.0;

	double qpr = re_q;
	double qpi = im_q;

	double qpmod = re_q*re_q+im_q*im_q;
	if (1.0 <= qpmod) return;

	for (n=1; n<max_terms; n++)
	{
		/* compute n^(-s) */
		double term = pow (n, -s);

		/* compute 1/(q^n-1) */
		double qmr = qpr-1.0;
		double qmi = qpi;
		double qm = qmr*qmr + qmi*qmi;

		double qvr = qmr/qm;
		double qvi = -qmi/qm;

		rep += qvr *term;
		imp += qvi *term;

		/* compute q^n */
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

	*prep = rep;
	*pimp = imp;
}

static double plouffe_series (double re_q, double im_q, int itermax, double param)
{
	ess = param;
	ess = 3.0;
	double rep, imp;
	plouffe_series_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	return rep;
}

DECL_MAKE_HISTO(plouffe_series);

/* --------------------------- END OF LIFE ------------------------- */

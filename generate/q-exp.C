/*
 * q-exp.C
 *
 * FUNCTION:
 * display q-exponential (degenerate basic hypergeoetric series)
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
#include "totient.h"


static int max_terms;

static void q_exp_c (double re_q, double im_q, double *prep, double
*pimp, double zee)
{
	int i;
	*prep = 0.0;
	*pimp = 0.0;

	double rep = 1.0;
	double imp = 0.0;

	/* power of q */
	double qpr = 1.0;
	double qpi = 0.0;

	double qpmod = re_q*re_q+im_q*im_q;
	if (1.0 <= qpmod) return;

	for (i=0; i<max_terms; i++)
	{
		/* compute 1-zq^k */
		double red = 1.0 - zee * qpr;
		double imd = -zee * qpi;

		/* compute 1/(1-zq^k) */
		double deno = 1.0/(red*red+imd*imd);
		double ret = red * deno;
		double imt = -imd * deno;

		/* compute prod * 1/(1-zq^k) */
		double tmp = rep*ret - imp * imt;
		imp = rep*imt + imp* ret;
		rep = tmp;

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

	*prep = rep;
	*pimp = imp;
}

static double q_exp_series (double re_q, double im_q, int itermax, double param)
{
	max_terms = itermax;
	double rep, imp;
	q_exp_c (re_q, im_q, &rep, &imp, param);
	// return sqrt (rep*rep+imp*imp);
	return rep;
}

DECL_MAKE_HISTO(q_exp_series);

/* --------------------------- END OF LIFE ------------------------- */

/*
 * zeta.C
 *
 * FUNCTION:
 * display binomial sums over the Riemann zeta
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

#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_psi.h>

#include "brat.h"

// ======================================================
// brute-force factorial function
inline long double factorial (int n)
{
	int k;
	long double fac = 1.0L;
	for (k=2; k<=n; k++)
	{
		fac *= (long double) k;
	}
	if (0>n) fac = 0.0L;
	return fac;
}

// ======================================================
// complex-valued binomial coefficent
// must have m>=0
// returns z*(z-1)*(z-2)...*(z-m+1) / m!
//
inline long double complex cbinomial (long double complex z, int m)
{
	if(0>m) return 0.0L;

	int k;
	long double complex bin = 1.0L;
	long double fac = 1.0L;
	for (k=1; k<=m; k++)
	{
		bin *= z;
		// printf ("bin term k=%d bin=(%Lg,%Lg) z=(%Lg,%Lg)\n", 
		//    k, creall(bin), cimagl(bin), creall(z), cimagl(z));
		z -= 1.0L;
		fac *= (long double) k;

		// avoid exponent overflows with periodic divisions
		if (1.0e300 < fac)
		{
			bin /= fac;
			fac = 1.0L;
		}
	}
	bin /= fac;
	return bin;
}
#endif

static void zeta_series_c (double re_q, double im_q, double *prep, double *pimp)
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
		double t=1;

		rep += qpr *t;
		imp += qpi *t;

		qpmod = qpr*qpr + qpi*qpi;
		if (qpmod < 1.0e-30) break;
	}

	*prep = rep;
	*pimp = imp;
}

static double zeta_series (double re_q, double im_q, int itermax, double param)
{
	max_terms = itermax;
	double rep, imp;
	zeta_series_c (re_q, im_q, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	return rep;
}

DECL_MAKE_HISTO(zeta_series);

/* --------------------------- END OF LIFE ------------------------- */

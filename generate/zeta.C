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
void cbinomial (double zre, double zim, int m, double *pbre, double *pbim)
{
	if(0>m) return;

	int k;
	double bre = 1.0;
	double bim = 0.0;
	long double fac = 1.0L;
	for (k=1; k<=m; k++)
	{
		double tmp = zre * bre - zim*bim;
		bim = zre * bim + zim * bre;
		bre = tmp;

		zre -= 1.0;
		fac *= (long double) k;
	}
	bre /= fac;
	bim /= fac;

	*pbre = bre;
	*pbim = bim;
}

static void zeta_series_c (double re_q, double im_q, double *prep, double *pimp)
{
	int i;
	*prep = 0.0;
	*pimp = 0.0;

	double rep = 0.0;
	double imp = 0.0;

	for (i=2; i<40; i++)
	{
		double bre, bim;
		cbinomial (re_q, im_q, i, &bre, &bim);

		double zetam1 = gsl_sf_zetam1 (i);
		if (i%2)
		{
			rep -= bre * zetam1;
			imp -= bim * zetam1;
		}
		else
		{
			rep += bre * zetam1;
			imp += bim * zetam1;
		}
	}

	*prep = rep;
	*pimp = imp;
}

static double zeta_series (double re_q, double im_q, int itermax, double param)
{
	double rep, imp;
	double re_z, im_z;

	double deno = re_q - 1.0;
	deno = deno*deno + im_q*im_q;
	deno = 1.0/deno;
	re_z = 2.0*im_q* deno;
	im_z = re_q*re_q + im_q*im_q - 1.0;
	im_z *= deno;

	zeta_series_c (re_z, im_z, &rep, &imp);
	// return sqrt (rep*rep+imp*imp);
	// return rep;
	// return imp;
	return (atan2 (imp,rep)+M_PI)/(2.0*M_PI);
}

DECL_MAKE_HISTO(zeta_series);

/* --------------------------- END OF LIFE ------------------------- */

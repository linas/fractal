
/*
 * gkw-integrand.c
 *
 * FUNCTION:
 * examine the integrand of the GKW integral to find the 
 * saddle point(s).
 *
 * Linas Vepstas January 2010
 */


#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "binomial.h"
#include "harmonic.h"
#include "zmp.h"


// Return the integrand 
void
grand(cpx_t prod, int m, int n, cpx_t s, int prec)
{
	cpx_t m2s, tms, sn1;

	cpx_init (m2s);
	cpx_init (tms);
	cpx_init (sn1);

	// sn1 = s+n+1
	cpx_add_ui (sn1, s, n+1, 0);

	// tms = 2-s
	cpx_neg (tms, s);
	cpx_add_ui (tms, tms, 2, 0);

	// m2s = m+2-s
	cpx_add_ui (m2s, tms, m, 0);

	cpx_t gm2s, gtms, gsn1, gs;
	cpx_init (gs);
	cpx_init (gm2s);
	cpx_init (gtms);
	cpx_init (gsn1);

	cpx_gamma(gs, s, prec);
	cpx_gamma(gm2s, m2s, prec);
	cpx_gamma(gtms, tms, prec);
	cpx_gamma(gsn1, sn1, prec);

	cpx_t zeta;
	cpx_init (zeta);
	cpx_borwein_zeta(zeta, m2s, prec);

	cpx_sub_ui(prod, zeta, 1, 0);
	cpx_mul(prod, prod, gm2s);
	cpx_mul(prod, prod, gs);
	cpx_div(prod, prod, gtms);
	cpx_div(prod, prod, gsn1);
}

long double
grand_modulus(int m, int n, cpx_t s, int prec)
{
	cpx_t prod;
	cpx_init (prod);

	grand(prod, m, n, s, prec);

	mpf_t modulus;
	mpf_init (modulus);
	cpx_abs(modulus, prod);

	return mpf_get_d (modulus);
}

long double
grand_phase(int m, int n, cpx_t s, int prec)
{
	cpx_t prod;
	cpx_init (prod);

	grand(prod, m, n, s, prec);

	mpf_t phase;
	mpf_init (phase);
	fp_arctan2 (phase, prod[0].im, prod[0].re, prec);

	double f = mpf_get_d (phase);
	f += M_PI;
	f /= 2.0*M_PI;
	return f;
}


static double gkw_integrand (double x, double y, int itermax, double param)
{
	static int init = 0;
	int prec = 60;
	if (0 == init)
	{
		/* Set the precision (number of binary bits) = prec*log(10)/log(2) */
		mpf_set_default_prec (3.3*prec);
		init = 1;
	}

	cpx_t s;
	cpx_init(s);
	cpx_set_d(s, x, y);

	int n = 10;
	int m = 10;
	// double gkw = grand_modulus(m,n,s, prec);
	double gkw = grand_phase(m,n,s, prec);

	printf("gr(%f +i %f) = %g\n", x, y, gkw);

	return gkw;
}

DECL_MAKE_HEIGHT(gkw_integrand);

/* --------------------------- END OF LIFE ------------------------- */

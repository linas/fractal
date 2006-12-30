/*
 * polylog.C
 *
 * High-precison Polylogarithm function, using the 
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas December 2006
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "../misc/bignum/zmp.h"

/* ============================================================================= */

static cpx_t zeta, ess, zee;
static int prec;

static void psi_init (void)
{
	/* the decimal precison (number of decimal places) */
	// prec = 90;
	// prec = 60;
	prec = 45;
	// prec = 25;

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
	int bits = (int) (v + 30);
	
	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	cpx_init (zeta);
	cpx_init (ess);
	cpx_init (zee);

	cpx_set_d (ess, 0.5, 14.134725141734);
	cpx_set_d (ess, 0.5, 21.022039638771);
	cpx_set_d (ess, 0.5, 37.5861781588256712);
}
	
static double plogger (double re_q, double im_q, int itermax, double param)
{
	static int init = 0;
	if (!init) {psi_init(); init=1; }
		  
	// printf ("duude compute %g  %g \n", re_z, im_z);

	cpx_set_d (zee, re_q, im_q);
	cpx_polylog (zeta, ess, zee, prec);

	double frea = mpf_get_d (zeta[0].re);
	double fima = mpf_get_d (zeta[0].im);

	double fmag = sqrt(frea*frea + fima*fima);
	return fmag;

	double phase = atan2 (fima, frea);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

DECL_MAKE_HEIGHT(plogger);

/* --------------------------- END OF LIFE ------------------------- */

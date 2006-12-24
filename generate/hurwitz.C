/*
 * hurwitz.C
 *
 * High-precison Hurwitz zeta function, using the 
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

static cpx_t zeta, ess;
mpf_t que;
static int prec;

static void psi_init (void)
{
	/* the decimal precison (number of decimal places) */
	prec = 50;

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
	/* XXX this should be stirling's approx for binomial */ 
	int bits = (int) (v + 30);
	
	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	mpf_init (que);
	cpx_init (zeta);
	cpx_init (ess);

	mpf_set_d (ess[0].re, 0.5);
}
	
static double hurl (double re_q, double im_q, int itermax, double param)
{
	static int init = 0;
	if (!init) {psi_init(); init=1; }
		  
	// printf ("duude compute %g  %g \n", re_z, im_z);
	mpf_set_d (que, re_q);
	mpf_set_d (ess[0].im, im_q);

	cpx_periodic_zeta (zeta, ess, que, prec);

	double frea = mpf_get_d (zeta[0].re);
	double fima = mpf_get_d (zeta[0].im);

	double fmag = sqrt(frea*frea + fima*fima);
	return fmag;
	
	double phase = atan2 (fima, frea);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

DECL_MAKE_HEIGHT(hurl);

/* --------------------------- END OF LIFE ------------------------- */

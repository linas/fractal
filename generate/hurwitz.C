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
	// prec = 90;
	// prec = 60;
	// prec = 40;
	prec = 15;

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
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

	// Its not sqrt, sqrt is too strong
	// im_q *= sqrt (re_q);
	// im_q *= pow (re_q, 1.0/3.0);
	// im_q *= pow (re_q, re_q); // nasty ugly turns down hard
	// im_q *= sqrt(sqrt (re_q));  // -- this is too much
	// im_q *= sqrt(sqrt(sqrt (re_q)));  // -- this is not enough
	// im_q *= pow (re_q, 1.0/6.0);  // not enough
	// im_q *= pow (re_q, 1.0/5.0);  // not enough  -- if turnup, then not obvious at lo res
	// im_q *= pow (re_q, 1.0/4.0);  // not enough .. but getting there -- subtle upturns at end
	// im_q *= pow (re_q, 1.0/3.0);  // Hmm. not enough in midrange, but upturms
	// im_q *= pow (re_q, 1.0/3.333333333);  // upturns
	// 
	// im_q *= pow (re_q, 24.0/(50.0*im_q)); // still sloping up, but not to hard
	// im_q *= pow (re_q, 12.0/(50.0*im_q));  // not strong enough, but still upturns
	// im_q *= pow (re_q, 2.0/sqrt(50.0*im_q));  //strong enough, upturning
	// im_q *= pow (re_q, 1.5/log(50.0*im_q));  // ohhh almost flat, but still upturning 
	// im_q *= pow (re_q, 1.2/log(50.0*im_q));  // down up,
	im_q *= pow (re_q, 1.2/log(100.0*im_q));  // down up,
	// im_q *= pow (re_q, 1.0/log(50.0*im_q));     // down 

	mpf_set_d (ess[0].im, 100.0*im_q);

	cpx_periodic_zeta (zeta, ess, que, prec);
	// cpx_hurwitz_zeta (zeta, ess, que, prec);

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

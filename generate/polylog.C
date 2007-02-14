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

static cpx_t zeta, ess, zee,z2, ph;
static int prec;

static void psi_init (int cmd_prec, double ims)
{
	/* the decimal precison (number of decimal places) */
	prec = 90;
	// prec = 60;
	// prec = 55;
	// prec = 45;
	// prec = 35;
	// prec = 25;
	prec = 20;
	// prec= ims;

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);
	
	/* the variable-precision calculations are touchy about this */
	int bits = (int) (v + 30 + 3.14*ims);
	
	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	cpx_init (zeta);
	cpx_init (ess);
	cpx_init (zee);
	cpx_init (z2);
	cpx_init (ph);

	// cpx_set_d (ess, 0.5, 14.134725141734);
	// cpx_set_d (ess, 0.5, 21.022039638771);
	// cpx_set_d (ess, 0.5, 37.5861781588256712);
	// cpx_set_d (ess, 0.5, 240.0);
	cpx_set_d (ess, 0.5, ims);

	cpx_set_d (zee, 0.0, 0.1);
	cpx_polylog_sheet_g0_action (ph, ess, 1, prec);
}

static double plogger (double re_q, double im_q, int itermax, double param)
{
	static int init = 0;
	if (!init) {psi_init(itermax, param); init=1; }
		  
	//printf ("duude compute %g  %g \n", re_q, im_q);

#ifdef S_PLANE
	cpx_set_d (ess, re_q, im_q);
	cpx_polylog (zeta, ess, zee, prec);
#endif

#define TRIPLICATION
#ifdef TRIPLICATION
	cpx_set_d (zee, re_q, im_q);
	cpx_polylog (zeta, ess, zee, prec);

	cpx_set_d (ph, -0.5, 0.866025404);
	cpx_mul (zee,zee,ph);
	cpx_polylog (z2, ess, zee, prec);
	cpx_add (zeta, zeta, z2);

	cpx_mul (zee,zee,ph);
	cpx_polylog (z2, ess, zee, prec);
	cpx_add (zeta, zeta, z2);
#endif

#ifdef DIRICHLET_CHARACTER_3_2
	cpx_set_d (zee, re_q, im_q);
	cpx_set_d (ph, -0.5, 0.866025404);
	cpx_mul (zee,zee,ph);
	cpx_polylog (zeta, ess, zee, prec);

	cpx_mul (zee,zee,ph);
	cpx_polylog (z2, ess, zee, prec);
	cpx_sub (zeta, zeta, z2);
#endif

#ifdef DIRICHLET_CHARACTER_4_2
	cpx_set_d (zee, re_q, im_q);
	cpx_times_i (zee,zee);
	cpx_polylog (zeta, ess, zee, prec);

	cpx_neg (zee,zee);
	cpx_polylog (z2, ess, zee, prec);
	cpx_sub (zeta, zeta, z2);
#endif

#ifdef DIRICHLET_CHARACTER_5_2
	cpx_set_d (zee, re_q, im_q);
	cpx_set_d (ph, cos(M_PI*0.4), sin (M_PI*0.4));
	cpx_mul (zee,zee,ph);
	cpx_polylog (zeta, ess, zee, prec);

	cpx_mul (zee,zee,ph);
	cpx_polylog (z2, ess, zee, prec);
	cpx_sub (zeta, zeta, z2);

	cpx_mul (zee,zee,ph);
	cpx_polylog (z2, ess, zee, prec);
	cpx_sub (zeta, zeta, z2);

	cpx_mul (zee,zee,ph);
	cpx_polylog (z2, ess, zee, prec);
	cpx_add (zeta, zeta, z2);
#endif

#ifdef SAME_AS_THREE_NEG_G1S
	cpx_polylog_sheet (z2, ess, zee, 0, -3, prec);
	cpx_sub (zeta, zeta, z2);
#endif

#ifdef THREE_NEG_G1S
	cpx_polylog_sheet_g1_action (z2, ess, zee, 0, -1, prec);
	cpx_sub (zeta, zeta, z2);
	cpx_polylog_sheet_g1_action (z2, ess, zee, -1, -1, prec);
	cpx_sub (zeta, zeta, z2);
	cpx_polylog_sheet_g1_action (z2, ess, zee, -2, -1, prec);
	cpx_sub (zeta, zeta, z2);
#endif

#if OTHER
	cpx_polylog_sheet_g1_action(z2, ess, zee, 0, 1, prec);
	cpx_mul (z2,z2,ph);
	cpx_sub(zeta, zeta, z2);
	cpx_polylog_sheet_g1_action(z2, ess, zee, 0, -1, prec);
	cpx_sub(zeta, zeta, z2);
	// cpx_polylog_sheet_g1_action(z2, ess, zee, -1, -1, prec);
	// cpx_sub(zeta, zeta, z2);
#endif

	double frea = mpf_get_d (zeta[0].re);
	double fima = mpf_get_d (zeta[0].im);

	// double fmag = sqrt(frea*frea + fima*fima);
	//	return fmag;

	double phase = atan2 (fima, frea);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

DECL_MAKE_HEIGHT(plogger);

/* --------------------------- END OF LIFE ------------------------- */

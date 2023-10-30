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
#include "zmp.h"

/* =================================================================== */

static cpx_t zeta, ess, zee, z2, ph;
static int prec;
static int branch = 0;

// Additional paramters are:
// sigma (float, real part of s)
// tau (float, imaginary part of s)
// branch (integer, 0 = main, see code for details)
// prec (integer, decimal precision)

static void psi_init ()
{
	// Sigma passed as parameter
	double sigma = 0.5;
	if (0 < param_argc)
		sigma = atof(param_argv[0]);
	double tau = 14;
	if (1 < param_argc)
		tau = atof(param_argv[1]);

	// branch to explore
	if (2 < param_argc)
		branch = atoi(param_argv[2]);

	/* The decimal precison (number of decimal places).
	 * prec=20 generates nice images. Some caution needed;
	 * if set too low, artifacts show up.   For example,
	 * setting prec=20 causes visible artifacts on the g1
	 * sheet at s = 0.5 +i 27 and above.
	 * The correct solution seems to be to add 1/3 of a unit
	 * of decimal precision for each additional bump in tau.
	 * The routines seem quite touchy about this.
	 * Hand-tuned; seems that this gives OK performance w/o
	 * loss of picture quality.
	 *
	 * prec=25 causes causes visible artifacts on the g1
	 * sheet at s = 0.5 +i 68 and above.
	 */
	prec = 25 + (int) (0.33 * tau);
	prec = 25 + (int) (0.5 * tau);
	prec = 35 + (int) (0.5 * tau);

	// Precision passed as parameter. Over-rides default guess above.
	if (3 < param_argc)
		prec = atoi(param_argv[3]);

	/* Compute number of binary bits this corresponds to. */
	/* log(10.0) / log(2.0) == 3.322 */
	double v = 3.322 * ((double) prec);
	int bits = (int) (v + 30 + 3.322*tau);

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

	cpx_set_d (ess, sigma, tau);

	printf("sigma=%g tau=%g branch=%d prec=%d bits=%d\n",
		sigma, tau, branch, prec, bits);

	cpx_set_d (zee, 0.0, 0.1);

	// cpx_polylog_sheet_g0_action returns -exp(-2pi is)
	// which we use multiplicatively to go around the branch at zero
	cpx_polylog_sheet_g0_action (ph, ess, 1, prec);
}

static double plogger (double re_q, double im_q, int itermax, double param)
{
	static int init = 0;
	if (!init) { psi_init(); init=1; }

	// printf ("duude compute %g  %g\n", re_q, im_q);

#ifdef S_PLANE
	cpx_set_d (ess, re_q, im_q);
	cpx_polylog (zeta, ess, zee, prec);
#endif

#ifdef TRIPLICATION
	cpx_set_d (zee, re_q, im_q);
	cpx_polylog (zeta, ess, zee, prec);

	// sqrt(3) = 0.866025404 so this is at 120 degrees
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

#define ZEE_PLANE
#ifdef ZEE_PLANE
	// This is the ordinary, generic polylog on the z plane, used
	// for the critical line movie.
	cpx_set_d (zee, re_q, im_q);
	cpx_polylog (zeta, ess, zee, prec);

	// Wind around the z=+1 cut in the left-hand direction
	if (-1 == branch)
	{
		cpx_polylog_sheet_g1_action(z2, ess, zee, 0, -1, prec);
		cpx_sub(zeta, zeta, z2);
	}

	// Wind around the z=+1 cut in the right-hand direction
	if (1 == branch)
	{
		cpx_polylog_sheet_g1_action(z2, ess, zee, 0, 1, prec);
		cpx_sub(zeta, zeta, z2);
	}

	// Wind around the z=+1 cut in the left-hand direction
	// Next, around the z=0 cut in the right-hand direction.
	// This is janky, but seems to be correct.
	if (2 == branch)
	{
		mpf_neg(zee[0].im, zee[0].im);
		cpx_polylog_sheet_g1_action(z2, ess, zee, 0, -2, prec);
		cpx_sub(zeta, zeta, z2);
	}

	// Wind around the z=+1 cut in the left-hand direction
	// Place the cut from z=0 so it's aimed to the left, not
	// the right. This gives a nice view of the double-cut.
	if (3 == branch)
	{
		if (mpf_cmp_ui(zee[0].im, 0) <= 0)
		{
			cpx_polylog_sheet_g1_action(z2, ess, zee, 0, -1, prec);
		}
		else
		{
			mpf_neg(zee[0].im, zee[0].im);
			cpx_polylog_sheet_g1_action(z2, ess, zee, 0, -2, prec);
		}
		cpx_sub(zeta, zeta, z2);
	}
#endif

	double frea = cpx_get_re(zeta);
	double fima = cpx_get_im(zeta);

	// double fmag = sqrt(frea*frea + fima*fima);
	//	return fmag;

	double phase = atan2 (fima, frea);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

DECL_MAKE_HEIGHT(plogger);

/* --------------------------- END OF LIFE ------------------------- */

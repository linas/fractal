/*
 * gamma.C
 *
 * High-precison Gamma function, using the 
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas December 2006
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "coord-xforms.h"
#include "../misc/bignum/zmp.h"

/* ======================================================================= */

static cpx_t gam, z;
static int prec;

static void psi_init (void)
{
	/* the decimal precison (number of decimal places) */
	prec = 120;

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
	/* XXX this should be stirling's approx for binomial */ 
	int bits = (int) (v + 30);
	
	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	cpx_init (gam);
	cpx_init (z);
}
	
static double cgamma (double re_q, double im_q, int itermax, double param)
{
	static int init = 0;
	if (!init) {psi_init(); init=1; }
		  
#if 0
	// double mag = 1.0 - sqrt (re_q*re_q + im_q*im_q);
	double mag = 1.0 - (re_q*re_q + im_q*im_q);

	if (mag <= 0.03) return 0.0;

	double re_z = re_q / mag;
	double im_z = im_q / mag;
#endif

	double re_z = re_q;
	double im_z = im_q;

#if USE_POINCARE_DISK
	double mag = (re_q*re_q + im_q*im_q);
	if (mag >= 0.97) return 0.0;
	poincare_disk_to_plane_coords (re_q, im_q, &re_z, &im_z);
#endif

	cpx_set_d (z, re_z, im_z);

	// printf ("duude compute %g  %g \n", re_z, im_z);
	cpx_gamma (gam, z, prec);

	double frea = mpf_get_d (gam[0].re);
	double fima = mpf_get_d (gam[0].im);

	double phase = atan2 (fima, frea);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

DECL_MAKE_HEIGHT(cgamma);

/* --------------------------- END OF LIFE ------------------------- */

/*
 * sho.c
 *
 * High-precison Simple Harmonic Oscillator, using the 
 * Gnu Multiple-precision library.
 *
 * Actually, high-precision eigenfunctions for arbitrary eigenvalue
 * 
 * Linas Vepstas November 2006
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "coord-xforms.h"
#include "../misc/bignum/zmp.h"

/* ======================================================================= */
/* compute psi_one eigenfunction for complex-valued s
 */

void psi_one (mpf_t re_psi, mpf_t im_psi, 
              double re_y, double im_y, unsigned int prec, double param)
{
	double relam = 1.5;
	double imlam = 0.0;

	relam = param;
	
	cpx_t em, a, b, z, ex;
	cpx_init (em);
	cpx_init (a);
	cpx_init (b);
	cpx_init (z);
	cpx_init (ex);

	/* a = 1/4 - lambda/2 */
	cpx_set_d (a, 0.25-0.5*relam, -0.5*imlam);

	/* b=1/2 */
	cpx_set_ui (b, 1, 0);
	cpx_div_ui (b, b, 2);
	
	/* z = y^2 */
	double reysq = re_y*re_y - im_y*im_y;
	double imysq = 2.0*re_y*im_y;
	cpx_set_d(z, reysq, imysq);
	
	cpx_confluent (em, a, b, z, prec);

	/* exp (-y^2/2) */
	cpx_neg (z, z);
	cpx_div_ui (z, z, 2);
	cpx_exp (ex, z, prec);

	/* psi = exp (-y^2/2) * M(a,b,y^2) */
	cpx_mul (em, em, ex);

	mpf_set (re_psi, em[0].re);
	mpf_set (im_psi, em[0].im);
						 
	cpx_clear (em);
	cpx_clear (a);
	cpx_clear (b);
	cpx_clear (z);
	cpx_clear (ex);

}

/* ============================================================================= */

static mpf_t re_a, im_a;
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
	
	mpf_init (re_a);
	mpf_init (im_a);
}
	
static double psi (double re_q, double im_q, int itermax, double param)
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

	// printf ("duude compute %g  %g \n", re_z, im_z);
	psi_one (re_a, im_a, re_z, im_z, prec, param);

	double frea = mpf_get_d (re_a);
	double fima = mpf_get_d (im_a);

	double fmag = sqrt(frea*frea + fima*fima);
	return fmag;
	
	double phase = atan2 (fima, frea);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

DECL_MAKE_HEIGHT(psi);

/* --------------------------- END OF LIFE ------------------------- */

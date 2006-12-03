
/*
 * Simple Harmonic Oscilltor
 *
 * Linas November 2006
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include "mp-complex.h"
#include "mp-hyper.h"
#include "mp-misc.h"
#include "mp-trig.h"

/* ======================================================== */
/* The two eigenfunctions of the simple harmonic oscillator */

void psi_1 (cpx_t psi, cpx_t lambda, cpx_t y, int prec)
{
	cpx_t a, b, z;
	cpx_init (a);
	cpx_init (b);
	cpx_init (z);

	/* b=1/2 */
	cpx_set_ui (b, 1, 0);
	mpf_div_ui (b[0].re, b[0].re, 2);

	/* a= 1/4-lambda/2 */
	cpx_set (a,b);
	cpx_sub (a, a, lambda);
	cpx_div_ui (a, a, 2);

	/* z = y^2 */
	cpx_mul (z, y, y);

	cpx_confluent (psi, a, b, z, prec);

	/* psi_1 = exp (-y^2/2) * M(a,b,z) */
	cpx_div_ui (z, z, 2);
	cpx_neg (z, z);
	cpx_exp (z, z, prec);
	cpx_mul (psi, psi, z);
	
	cpx_clear (a);
	cpx_clear (b);
	cpx_clear (z);
}

void psi_2 (cpx_t psi, cpx_t lambda, cpx_t y, int prec)
{
	cpx_t a, b, z;
	cpx_init (a);
	cpx_init (b);
	cpx_init (z);

	/* b=3/2 */
	cpx_set_ui (b, 3, 0);
	mpf_div_ui (b[0].re, b[0].re, 2);

	/* a= (3/2-lambda)/2 */
	cpx_set (a,b);
	cpx_sub (a, a, lambda);
	cpx_div_ui (a, a, 2);

	/* z = y^2 */
	cpx_mul (z, y, y);

	cpx_confluent (psi, a, b, z, prec);

	/* psi_2 = y * exp (-y^2/2) * M(a,b,z) */
	cpx_div_ui (z, z, 2);
	cpx_neg (z, z);
	cpx_exp (z, z, prec);
	cpx_mul (psi, psi, z);
	cpx_mul (psi, psi, y);
	
	cpx_clear (a);
	cpx_clear (b);
	cpx_clear (z);
}

/* ======================================================== */
/* cheap hack for the gamma function. Incorrect, but holds water.
 * Approximates gamma(z)=1 for 1<z<2.
 */
static void gamma_hack (cpx_t gam, cpx_t z)
{
	double flo = mpf_get_d (z[0].re);
	unsigned int intpart = (unsigned int) floor (flo);
	intpart ++;
	cpx_set (gam, z);
	mpf_sub_ui (gam[0].re, gam[0].re, inpart);
	cpx_poch_rising (gam, gam, intpart);
}

void eta_1 (cpx_t psi, cpx_t lambda, cpx_t y, int prec)
{
	cpx_t psi;
	cpx_t pha;
	cpx_init (psi);
	cpx_init (pha);

	psi_1 (psi, lambda, y, prec);

	mpf_t pi;
	mpf_init (pi);
	fp_pi (pi, prec);

	/* phase term, exp (i pi lambda/2) */
	cpx_set (pha, lambda);
	cpx_div_ui (pha, 2);
	cpx_mul_mpf (pha, pi);
	cpx_times_i (pha, pha);
	cpx_exp (pha, pha);
	cpx_mul (psi, psi, pha);

	/* power of two term */
	mpf_set_ui (pi, 2);
	cpx_set (pha, lambda);
	cpx_div_ui (pha, 2);
	cpx_pow_mpf (pha, pi, pha, prec);
	

	cpx_clear (psi);
	cpx_clear (pha);
	mpf_clear (pi);
}

/* ======================================================== */
int
main (int argc, char *argv[])
{
	int prec = 250;
	
	double lam;
	lam = atof (argv[1]);

	cpx_t ps1, ps2, lambda, z;
	cpx_init (ps1);
	cpx_init (ps2);
	cpx_init (lambda);
	cpx_init (z);

	cpx_set_ui (lambda, 1,0);
	mpf_set_d (lambda[0].re, lam);

	cpx_t ex;
	cpx_init (ex);
	cpx_set_ui (ex, 0,0);

	double x;
	for (x=-6.0; x<6.01; x+=0.1)
	{
		mpf_set_d (ex[0].re, x);

		psi_1 (ps1, lambda, ex, prec);
		psi_2 (ps2, lambda, ex, prec);

		printf ("%g", x);
		fp_prt ("\t", ps1[0].re);
		fp_prt ("\t", ps2[0].re);
		printf ("\n");
	}

	return 0;
}

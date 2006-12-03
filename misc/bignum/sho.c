
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

void psi_1 (cpx_t *psi, cpx_t *lambda, cpx_t *y, int prec)
{
	cpx_t a, b, z;
	cpx_init (&a);
	cpx_init (&b);
	cpx_init (&z);

	/* b=1/2 */
	cpx_set_ui (&b, 1, 0);
	mpf_div_ui (b.re, b.re, 2);

	/* a= 1/4-lambda/2 */
	cpx_set (&a,&b);
	cpx_sub (&a, &a, lambda);
	cpx_div_ui (&a, &a, 2);

	/* z = y^2 */
	cpx_mul (&z, y, y);

	cpx_confluent (psi, &a, &b, &z, prec);

	/* psi_1 = exp (-y^2/2) * M(a,b,z) */
	cpx_div_ui (&z, &z, 2);
	cpx_neg (&z, &z);
	cpx_exp (&z, &z, prec);
	cpx_mul (psi, psi, &z);
	
	cpx_clear (&a);
	cpx_clear (&b);
	cpx_clear (&z);
}

void psi_2 (cpx_t *psi, cpx_t *lambda, mpf_t x, int prec)
{
	cpx_t a, b, z;
	cpx_init (&a);
	cpx_init (&b);
	cpx_init (&z);

	/* b=3/2 */
	cpx_set_ui (&b, 3, 0);
	mpf_div_ui (b.re, b.re, 2);

	/* a= (3/2-lambda)/2 */
	cpx_set (&a,&b);
	cpx_sub (&a, &a, lambda);
	cpx_div_ui (&a, &a, 2);

	/* z = y^2 */
	mpf_set (z.re, x);
	mpf_mul (z.re, z.re, x);
	mpf_set_ui (z.im, 0);

	cpx_confluent (psi, &a, &b, &z, prec);

	/* psi_2 = y * exp (-y^2/2) * M(a,b,z) */
	mpf_div_ui (z.re, z.re, 2);
	mpf_neg (z.re, z.re);
	fp_exp (z.im, z.re, prec);
	cpx_mul_mpf (psi, psi, z.im);
	cpx_mul_mpf (psi, psi, x);
	
	cpx_clear (&a);
	cpx_clear (&b);
	cpx_clear (&z);
}

int
main (int argc, char *argv[])
{
	int prec = 250;
	
	double lam;
	lam = atof (argv[1]);

	cpx_t ps1, ps2, lambda, z;
	cpx_init (&ps1);
	cpx_init (&ps2);
	cpx_init (&lambda);
	cpx_init (&z);

	cpx_set_ui (&lambda, 1,0);
	mpf_set_d (lambda.re, lam);

	mpf_t ex;
	mpf_init (ex);

	double x;
	for (x=-6.0; x<6.01; x+=0.1)
	{
		mpf_set_d (ex, x);

		psi_1 (&ps1, &lambda, ex, prec);
		psi_2 (&ps2, &lambda, ex, prec);

		printf ("%g", x);
		fp_prt ("\t", ps1.re);
		fp_prt ("\t", ps2.re);
		printf ("\n");
	}

	return 0;
}

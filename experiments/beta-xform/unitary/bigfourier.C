/*
 * bigfourier.C
 *
 * Compute Fourier transform of convergents
 * Feb 2024
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"
#include "bigseries.c"

#include <anant/mp-trig.h>

// -------------------------------------------------------

#define NBITS 1200
char bitseq[NBITS];
bool initialized = false;

static void do_init(void)
{
	int bprec = 1400;
	mpf_set_default_prec(bprec);
	printf("#\n# Default prec=%d bits\n#\n", bprec);

	// Set beta to exactly 1.6
	mpf_t beta;
	mpf_init(beta);
	mpf_set_ui(beta, 16);
	mpf_div_ui(beta, beta, 10);

	gen_bitseq(beta, bitseq, NBITS);
	printf("#\n# ");
	for (int i=0; i<70; i++)
		printf("%d", bitseq[i]);
	printf("\n#\n");
}

static double fourier(double re_q, double im_q, int itermax, double param)
{
	if (false == initialized)
	{
		initialized = true;
		do_init();
	}

#if 0
	// Frequency on the vertical axis, order on the horizontal.
	double x = im_q;
	int order = 800 * re_q;
#endif

	// Frequency angular direction, order on the radial.
	double phi = atan2(im_q, re_q) / (2.0*M_PI);
	double x = phi;
	int order = 800 * sqrt(re_q*re_q + im_q*im_q);

	// zeta = exp (i 2pi x)
	cpx_t zeta;
	cpx_init(zeta);
	cpx_set_d(zeta, 0, 2*M_PI*x);
	cpx_exp(zeta, zeta, 50);

	cpx_t sum;
	cpx_init(sum);
	ebz(sum, zeta, bitseq, order);

	double re = cpx_get_re(sum);
	double im = cpx_get_im(sum);
	double mod = sqrt(re*re+im*im);
	mod /= sqrt(order);

	cpx_clear(zeta);
	cpx_clear(sum);

	return mod;
}

DECL_MAKE_HEIGHT(fourier);
   
/* --------------------------- END OF LIFE ------------------------- */

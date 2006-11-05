/*
 * mp-consts.c
 *
 * High-precison constants, using the 
 * Gnu Multiple-precision library.
 *
 * Also, high-precision values of the series a_n 
 * 
 * Linas Vepstas July 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include "mp-cache.h"
#include "mp-consts.h"
#include "mp-trig.h"

/* ======================================================================= */
/**
 * fp_pi - return pi=3.14159... 
 * @prec - number of decimal places of precision
 *
 * Uses simple, brute-force Machin formula
 */
void fp_pi (mpf_t pi, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_pi;

	if (precision >= prec)
	{
		mpf_set (pi, cached_pi);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_pi);
	}

	/* Simple-minded Machin formula */
	mpf_t tmp;
	mpf_init(tmp);
	mpf_set_ui (tmp, 1);
	mpf_div_ui (tmp, tmp, 5);
	fp_arctan (pi, tmp, prec);

	mpf_mul_ui (pi, pi, 4);

	mpf_set_ui (tmp, 1);
	mpf_div_ui (tmp, tmp, 239);
	fp_arctan (tmp, tmp, prec);

	mpf_sub (pi, pi, tmp);
	mpf_mul_ui (pi, pi, 4);
	mpf_clear (tmp);

	mpf_set (cached_pi, pi);
	precision = prec;
}

/**
 * fp_e_pi - return e^pi 
 * @prec - number of decimal places of precision
 *
 * Uses simple, low-brow formula
 */
void fp_e_pi (mpf_t e_pi, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_e_pi;

	if (precision >= prec)
	{
		mpf_set (e_pi, cached_e_pi);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_e_pi);
	}

	fp_pi (e_pi, prec);
	fp_exp (e_pi, e_pi, prec);

	mpf_set (cached_e_pi, e_pi);
	precision = prec;
}


/* fp_euler
 * return Euler-Mascheroni const
 */
static void fp_euler_mascheroni_compute (mpf_t gam, unsigned int prec)
{
	/* power value, goes as log log n */
	// double en = log (prec) / log (2.0);
	double en = log (prec);
	en = 1.442695041 * (en + log (en));
	int n = (int) en;
	
	mpf_t maxterm;
	mpf_init (maxterm);
	mpf_set_ui (maxterm, 1);

	mpf_t z_n, twon, term, tmp, fact;
	mpf_init (z_n);
	mpf_init (twon);
	mpf_init (term);
	mpf_init (tmp);
	mpf_init (fact);
	mpf_set_ui (twon, 1);
	mpf_mul_2exp(twon, twon, n);
	mpf_set (z_n, twon);
	mpf_set_ui (fact, 1);
	mpf_div_ui (fact, fact, 2);
	mpf_set_ui (gam, 1);

	unsigned int k=2;
	while(1)
	{
		fp_harmonic (tmp, k);
		mpf_mul (term, z_n, tmp);
		mpf_mul (term, term, fact);

		mpf_add (gam, gam, term);

		/* don't go no farther than this */
		if (mpf_cmp (term, maxterm) < 0) break;

		k ++;
		mpf_mul (z_n, z_n, twon);
		mpf_div_ui (fact, fact, k);
	}

	mpf_mul (gam, gam, twon);

	fp_exp (tmp, twon, prec);
	mpf_div (gam, gam, tmp);

	mpf_set_ui (tmp, 2);
	fp_log(tmp, tmp, prec);
	mpf_mul_ui (tmp, tmp, n);
	mpf_sub (gam, gam, tmp);
	
	mpf_clear (z_n);
	mpf_clear (twon);
	mpf_clear (term);
	mpf_clear (tmp);
	mpf_clear (fact);

	mpf_clear (maxterm);
}

void fp_euler_mascheroni (mpf_t gam, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_gam;

	if (precision >= prec)
	{
		mpf_set (gam, cached_gam);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_gam);
	}

	fp_euler_mascheroni_compute (gam, prec);
	mpf_set (cached_gam, gam);
	precision = prec;
}

/* =============================== END OF FILE =========================== */


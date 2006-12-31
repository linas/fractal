/*
 * mp-consts.c
 *
 * High-precison constants, using the 
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas July 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include "mp-binomial.h"
#include "mp-complex.h"
#include "mp-consts.h"
#include "mp-trig.h"
#include "mp-zeta.h"

/* ======================================================================= */
/**
 * fp_half_sqrt_three - return 0.5*sqrt(3)= 0.86602...
 */

void fp_half_sqrt_three (mpf_t sqt)
{
	static unsigned int init=0;
	static mpf_t cached_sqt;

	if (init)
	{
		mpf_set (sqt, cached_sqt);
		return;
	}
	mpf_init (cached_sqt);

	mpf_set_ui (sqt, 3);
	mpf_sqrt (sqt, sqt);
	mpf_div_ui (sqt, sqt, 2);
	mpf_set (cached_sqt, sqt);

	init =1;
}

/* ======================================================================= */
/**
 * fp_e - return e=2.718281828...
 * @prec - number of decimal places of precision
 *
 * Uses simple, brute-force summation
 */

extern void fp_exp_helper (mpf_t ex, const mpf_t z, unsigned int prec);

void fp_e (mpf_t e, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_e;

	if (precision >= prec)
	{
		mpf_set (e, cached_e);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_e);
	}

	mpf_t one;
	mpf_init (one);
	mpf_set_ui (one, 1);
	fp_exp_helper (cached_e, one, prec);
	mpf_set (e, cached_e);

	mpf_clear (one);
	precision = prec;
}

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

/* ======================================================================= */
/**
 * fp_two_pi - return 2pi = 2 * 3.14159... 
 * @prec - number of decimal places of precision
 *
 * The idea is that it caches the value to avoid recomputation
 */
void fp_two_pi (mpf_t two_pi, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_two_pi;

	if (precision >= prec)
	{
		mpf_set (two_pi, cached_two_pi);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_two_pi);
	}

	fp_pi (two_pi, prec);
	mpf_mul_ui (two_pi, two_pi, 2);
	mpf_set (cached_two_pi, two_pi);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_pi_half - return pi/2 = 0.5 * 3.14159... 
 * @prec - number of decimal places of precision
 *
 * The idea is that it caches the value to avoid recomputation
 */
void fp_pi_half (mpf_t pih, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_pih;

	if (precision >= prec)
	{
		mpf_set (pih, cached_pih);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_pih);
	}

	fp_pi (pih, prec);
	mpf_div_ui (pih, pih, 2);
	mpf_set (cached_pih, pih);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_sqrt_two_pi - return sqrt(2pi) = sqrt (2 * 3.14159...) 
 * @prec - number of decimal places of precision
 *
 * The idea is that it caches the value to avoid recomputation
 */
void fp_sqrt_two_pi (mpf_t sqtpi, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_sqtpi;

	if (precision >= prec)
	{
		mpf_set (sqtpi, cached_sqtpi);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_sqtpi);
	}

	fp_two_pi (sqtpi, prec);
	mpf_sqrt (sqtpi, sqtpi);
	mpf_set (cached_sqtpi, sqtpi);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_log_two_pi - return log(2pi) = log(2 * 3.14159...) 
 * @prec - number of decimal places of precision
 *
 * The idea is that it caches the value to avoid recomputation
 */
void fp_log_two_pi (mpf_t ltp, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_ltp;

	if (precision >= prec)
	{
		mpf_set (ltp, cached_ltp);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_ltp);
	}

	fp_two_pi (ltp, prec);
	fp_log (ltp, ltp, prec);
	mpf_set (cached_ltp, ltp);
	precision = prec;
}

/* ======================================================================= */
/**
 * fp_log2 - return log(2)=0.69...
 * @prec - number of decimal places of precision
 *
 * Uses simple, brute-force summation
 */

void fp_log2 (mpf_t l2, unsigned int prec)
{
	static unsigned int precision=0;
	static mpf_t cached_log2;

	if (precision >= prec)
	{
		mpf_set (l2, cached_log2);
		return;
	}

	if (0 == precision)
	{
		mpf_init (cached_log2);
	}

	mpf_t two;
	mpf_init (two);
	mpf_set_ui (two, 2);
	fp_log (cached_log2, two, prec);
	mpf_set (l2, cached_log2);

	mpf_clear (two);
	precision = prec;
}

/* ======================================================================= */
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


/* ======================================================================= */
/** 
 * fp_euler - return Euler-Mascheroni const
 * @prec - number of decimal places of precision
 *
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

/* ======================================================================= */
/** 
 * fp_zeta_half - return zeta (1/2)
 * @prec - number of decimal places of precision
 *
 */
static void fp_zeta_half_compute (mpf_t gam, unsigned int prec)
{
	cpx_t ess, zeta;
	cpx_init (ess);
	cpx_init (zeta);
	
	mpf_set_d (ess[0].re, 0.5);
	mpf_set_ui (ess[0].im, 0);

	cpx_borwein_zeta (zeta, ess, prec);
	mpf_set (gam, zeta[0].re);

	cpx_clear (ess);
	cpx_clear (zeta);
}

void fp_zeta_half (mpf_t gam, unsigned int prec)
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

	fp_zeta_half_compute (gam, prec);
	mpf_set (cached_gam, gam);
	precision = prec;
}

/* =============================== END OF FILE =========================== */


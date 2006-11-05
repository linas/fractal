/*
 * mp-trig.c
 *
 * High-precison Elementary functions, using the 
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
 * i_pow - raise n to the m power
 */

void i_pow (mpz_t p, unsigned int n, unsigned int m)
{
	DECLARE_I_CACHE (cache);
	cache.disabled = 1;

	if ((1 == n) || (0 == m))
	{
		mpz_set_ui (p, 1); 
		return;
	}

	int hit = i_triangle_cache_check (&cache, n+m, m);
	if (hit)
	{
		i_triangle_cache_fetch (&cache, p, n+m, m);
	}
	else
	{
		i_pow (p, n, m-1);
		mpz_mul_ui (p, p, n);
		i_triangle_cache_store (&cache, p, n+m, m);
	}
}

/**
 * fp_inv_pow - raise n to the -m power, where m must be positive. 
 */

void fp_inv_pow (mpf_t p, unsigned int n, unsigned int m)
{
	DECLARE_FP_CACHE (cache);
	if (1 == n)
	{
		mpf_set_ui (p, 1); 
		return;
	}

	int hit = fp_triangle_cache_check (&cache, n+m, m);
	if (hit)
	{
		fp_triangle_cache_fetch (&cache, p, n+m, m);
	}
	else
	{
		mpz_t ip;
		mpz_init (ip);
		i_pow (ip, n, m);
		mpf_set_z (p, ip);
		mpf_ui_div (p, 1, p);
		mpz_clear (ip);
		fp_triangle_cache_store (&cache, p, n+m, m, 1);
	}
}

/* ======================================================================= */
/**
 * fp_exp -  Floating point exponential
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

static void fp_exp_helper (mpf_t ex, mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, fact, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (fact);
	mpf_init (term);

	mpf_set_ui (ex, 1);
	mpf_set_ui (fact, 1);
	mpf_set (zee, z);
	mpf_set (z_n, zee);
	
	// double mex = ((double) prec) * log (10.0) / log(2.0);
	double mex = ((double) prec) * 3.321928095;
	unsigned int imax = (unsigned int) (mex +1.0);
	mpf_t maxterm, one;
	mpf_init (maxterm);
	mpf_init (one);
	mpf_set_ui (one, 1);
	mpf_div_2exp (maxterm, one, imax);

	unsigned int n=1;
	while(1)
	{
		mpf_div (term, z_n, fact);
		mpf_add (ex, ex, term);
		
		/* don't go no farther than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n++;
		mpf_mul (z_n, z_n, zee);
		mpf_mul_ui (fact, fact, n);
	}
	
	mpf_clear (zee);
	mpf_clear (z_n);
	mpf_clear (fact);
	mpf_clear (term);

	mpf_clear (one);
	mpf_clear (maxterm);
}

void fp_exp (mpf_t ex, mpf_t z, unsigned int prec)
{
	if (mpf_cmp_ui (z, 0) > 0)
	{
		fp_exp_helper (ex, z, prec);
	}
	else
	{
		mpf_t zee;
		mpf_init (zee);
		mpf_neg (zee, z);
		fp_exp_helper (ex, zee, prec);
		mpf_ui_div (ex, 1, ex);
		mpf_clear (zee);
	}
}

/* ======================================================================= */
/**
 * fp_log_m1 -  Floating point logarithm
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

void fp_log_m1 (mpf_t lg, mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (term);

	mpf_set (zee, z);
	mpf_mul (z_n, zee, zee);
	mpf_set (lg, zee);
	
	// double mex = ((double) prec) * log (10.0) / log(2.0);
	double mex = ((double) prec) * 3.321928095;
	unsigned int imax = (unsigned int) (mex +1.0);
	mpf_t maxterm, one;
	mpf_init (maxterm);
	mpf_init (one);
	mpf_set_ui (one, 1);
	mpf_div_2exp (maxterm, one, imax);

	unsigned int n=2;
	while(1)
	{
		mpf_div_ui (term, z_n, n);
		mpf_add (lg, lg, term);
		
		/* don't go no farther than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n ++;
		mpf_mul (z_n, z_n, zee);
	}
	
	mpf_clear (zee);
	mpf_clear (z_n);
	mpf_clear (term);

	mpf_clear (one);
	mpf_clear (maxterm);
}

void fp_log (mpf_t lg, mpf_t z, unsigned int prec)
{
	mpf_t zee;
	mpf_init (zee);
	if (mpf_cmp_d(z, 1.5) > 0)
	{
		mpf_ui_div (zee, 1, z);
		mpf_ui_sub (zee, 1, zee);
		fp_log_m1 (lg, zee, prec);
	}
	else
	{
		mpf_ui_sub (zee, 1, z);
		fp_log_m1 (lg, zee, prec);
		mpf_neg (lg, lg);
	}
	mpf_clear (zee);
}

/* ======================================================================= */
/**
 * fp_arctan -  Floating point arctangent
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

void fp_arctan (mpf_t atn, mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, zsq, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (zsq);
	mpf_init (term);

	mpf_set (zee, z);
	mpf_mul (zsq, zee, zee);
	mpf_mul (z_n, zee, zsq);
	mpf_set (atn, zee);
	
	// double mex = ((double) prec) * log (10.0) / log(2.0);
	double mex = ((double) prec) * 3.321928095;
	unsigned int imax = (unsigned int) (mex +1.0);
	mpf_t maxterm, one;
	mpf_init (maxterm);
	mpf_init (one);
	mpf_set_ui (one, 1);
	mpf_div_2exp (maxterm, one, imax);

	unsigned int n=1;
	while(1)
	{
		mpf_div_ui (term, z_n, 2*n+1);
		if (n%2)
		{
			mpf_sub (atn, atn, term);
		}
		else
		{
			mpf_add (atn, atn, term);
		}
		
		/* don't go no farther than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n ++;
		mpf_mul (z_n, z_n, zsq);
	}
	
	mpf_clear (zee);
	mpf_clear (z_n);
	mpf_clear (zsq);
	mpf_clear (term);

	mpf_clear (one);
	mpf_clear (maxterm);
}

/* =============================== END OF FILE =========================== */


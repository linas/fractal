/*
 * mp-trig.c
 *
 * High-precison Elementary functions, using the 
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas July 2005
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include "mp-binomial.h"
#include "mp-cache.h"
#include "mp-complex.h"
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

static void fp_exp_helper (mpf_t ex, const mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, fact, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (fact);
	mpf_init (term);

	/* Make copy of argument now! */
	mpf_set (zee, z);
	mpf_set (z_n, zee);
	
	mpf_set_ui (ex, 1);
	mpf_set_ui (fact, 1);

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

void fp_exp (mpf_t ex, const mpf_t z, unsigned int prec)
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
 * fp_sine -  Floating point sine function
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

void fp_sine (mpf_t si, const mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, fact, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (fact);
	mpf_init (term);

	/* Make copy of argument now! */
	mpf_set (zee, z);
	mpf_set (z_n, zee);
	mpf_set_ui (si, 0);
	mpf_set_ui (fact, 1);
	
	// double mex = ((double) prec) * log (10.0) / log(2.0);
	prec += 2;
	double mex = ((double) prec) * 3.321928095;
	unsigned int imax = (unsigned int) (mex +1.0);
	mpf_t maxterm, one;
	mpf_init (maxterm);
	mpf_init (one);
	mpf_set_ui (one, 1);
	mpf_div_2exp (maxterm, one, imax);

	unsigned int n=1;
	unsigned int s=0;
	while(1)
	{
		mpf_div (term, z_n, fact);

		if (0 == s%2)
		{
			mpf_add (si, si, term);
		}
		else
		{
			mpf_sub (si, si, term);
		}
		
		/* don't go no farther than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n++;
		mpf_mul (z_n, z_n, zee);
		mpf_mul_ui (fact, fact, n);
		
		n++;
		mpf_mul (z_n, z_n, zee);
		mpf_mul_ui (fact, fact, n);

		s++;
	}
	
	mpf_clear (zee);
	mpf_clear (z_n);
	mpf_clear (fact);
	mpf_clear (term);

	mpf_clear (one);
	mpf_clear (maxterm);
}

/* ======================================================================= */
/**
 * fp_cosine -  Floating point cosine function
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

void fp_cosine (mpf_t co, const mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, fact, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (fact);
	mpf_init (term);

	/* Make copy of argument now! */
	mpf_set (zee, z);
	mpf_mul (z_n, zee, zee);
	mpf_set_ui (co, 1);
	mpf_set_ui (fact, 2);
	
	// double mex = ((double) prec) * log (10.0) / log(2.0);
	prec +=2;
	double mex = ((double) prec) * 3.321928095;
	unsigned int imax = (unsigned int) (mex +1.0);
	mpf_t maxterm, one;
	mpf_init (maxterm);
	mpf_init (one);
	mpf_set_ui (one, 1);
	mpf_div_2exp (maxterm, one, imax);

	unsigned int n=2;
	unsigned int s=1;
	while(1)
	{
		mpf_div (term, z_n, fact);

		if (0 == s%2)
		{
			mpf_add (co, co, term);
		}
		else
		{
			mpf_sub (co, co, term);
		}
		
		/* Don't go no farther than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n++;
		mpf_mul (z_n, z_n, zee);
		mpf_mul_ui (fact, fact, n);
		
		n++;
		mpf_mul (z_n, z_n, zee);
		mpf_mul_ui (fact, fact, n);

		s++;
	}
	
	mpf_clear (zee);
	mpf_clear (z_n);
	mpf_clear (fact);
	mpf_clear (term);

	mpf_clear (one);
	mpf_clear (maxterm);
}

/* ======================================================================= */
/**
 * cpx_exp -  complex-valued exp, built out of the real-valued funcs.
 */

void cpx_exp (cpx_t ex, const cpx_t const z, unsigned int prec)
{
	mpf_t mag, si, co;

	mpf_init (mag);
	mpf_init (si);
	mpf_init (co);

	fp_exp (mag, z->re, prec);
	fp_cosine (co, z->im, prec);
	fp_sine (si, z->im, prec);

	mpf_mul (ex->re, mag, co);
	mpf_mul (ex->im, mag, si);
	
	mpf_clear (mag);
	mpf_clear (si);
	mpf_clear (co);
}

/* ======================================================================= */
/**
 * fp_log_m1 -  Floating point logarithm
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

void fp_log_m1 (mpf_t lg, const mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (term);

	/* Make copy of argument now! */
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

void fp_log (mpf_t lg, const mpf_t z, unsigned int prec)
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

void fp_arctan (mpf_t atn, const mpf_t z, unsigned int prec)
{
	mpf_t zee, z_n, zsq, term;

	mpf_init (zee);
	mpf_init (z_n);
	mpf_init (zsq);
	mpf_init (term);

	/* Make copy of argument now! */
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

/* ======================================================================= */
/**
 * cpx_pow_mpf-- return q^s for complex s, real q.
 *
 * Brute-force algo, this thing is pretty slow, as it requires
 * a logarithm, an exp, sin and cos to be computed, each of which
 * are kinda slow ... 
 */
void cpx_pow_mpf (cpx_t powc, const mpf_t kq, const cpx_t const ess, int prec)
{
	mpf_t logkq, mag, pha;
	mpf_init (logkq);
	mpf_init (mag);
	mpf_init (pha);

	fp_log (logkq, kq, prec);
	
	/* magnitude is exp(re(s) * log(kq)) */
	mpf_mul (mag, ess->re, logkq);
	
	fp_exp (mag, mag, prec);

	/* phase is im(s) * log(kq)) */
	mpf_mul (pha, ess->im, logkq);

	fp_cosine (powc->re, pha, prec);
	mpf_mul (powc->re, mag, powc->re);
	
	fp_sine (powc->im, pha, prec);
	mpf_mul (powc->im, mag, powc->im);
	
	mpf_clear(logkq);
	mpf_clear(mag);
	mpf_clear(pha);
}

/* ======================================================================= */
/**
 * fp_pow_rc-- return (k+q)^s for complex s, integer k, real q.
 *
 * If q is held fixed, and k varied, then the values are cached,
 * allowing improved algorithm speeds.
 *
 * Overall, though, this thing is pretty slow, as it requires
 * a logarithm, an exp, sin and cos to be computed, each of which
 * are kinda slow ... 
 */
void fp_pow_rc (cpx_t powc, int k, const mpf_t q, const cpx_t const ess, int prec)
{
	DECLARE_FP_CACHE (re_powc);
	DECLARE_FP_CACHE (im_powc);
	static mpf_t cache_q;
	static int init = 0;

	if (!init)
	{
		init = 1;
		mpf_init (cache_q);
	}

	if (!mpf_eq(q,cache_q, prec*3.322))
	{
		fp_one_d_cache_clear (&re_powc);
		fp_one_d_cache_clear (&im_powc);
		mpf_set(cache_q,q);
	}

	if (prec <= fp_one_d_cache_check (&re_powc, k))
	{
		fp_one_d_cache_fetch (&re_powc, powc->re, k);
		fp_one_d_cache_fetch (&im_powc, powc->im, k);
		return;
	}
	
	mpf_t kq;
	mpf_init (kq);
	mpf_add_ui (kq, q, k);
	cpx_pow_mpf (powc, kq, ess, prec);
	mpf_clear (kq);

	fp_one_d_cache_check (&im_powc, k);
	fp_one_d_cache_store (&re_powc, powc[0].re, k, prec);
	fp_one_d_cache_store (&im_powc, powc[0].im, k, prec);
}

/* =============================== END OF FILE =========================== */


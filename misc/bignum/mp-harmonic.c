/* 
 * mp-harmonic.c
 *
 * Compute assorted harmonic functions to arbitrary precision
 *
 * Linas Vepstas October 2006
 */

#include <math.h>
#include <stdio.h>

#include <gmp.h>
#include "mp-cache.h"
#include "mp-complex.h"
#include "mp-harmonic.h"
#include "mp-trig.h"

/* =========================================================== */

/**
 * mp_pow_rc-- return (k+q)^s for complex s, integer k, real q.
 *
 * If q is held fixed, and k varied, then the values are cached,
 * allowing improved algorithm speeds.
 */
static void mp_pow_rc_helper (cpx_t *diri, int k, mpf_t q, cpx_t *ess, int prec)
{
	mpf_t kq, logkq, mag, pha;
	mpf_init (kq);
	mpf_init (logkq);
	mpf_init (mag);
	mpf_init (pha);

	mpf_add_ui (kq, q, k);
	
	fp_log (logkq, kq, prec);
	
	/* magnitude is exp(re(s) * log(k+q)) */
	mpf_mul (mag, ess->re, logkq);
	fp_exp (mag, mag, prec);

	/* phase is im(s) * log(k+q)) */
	mpf_mul (pha, ess->im, logkq);

	fp_cosine (diri->re, pha, prec);
	mpf_mul (diri->re, mag, diri->re);
	
	fp_sine (diri->im, pha, prec);
	mpf_mul (diri->im, mag, diri->im);
	
	mpf_clear(kq);
	mpf_clear(logkq);
	mpf_clear(mag);
	mpf_clear(pha);
}

void mp_pow_rc (cpx_t *diri, int k, mpf_t q, cpx_t *ess, int prec)
{
	DECLARE_FP_CACHE (re_diri);
	DECLARE_FP_CACHE (im_diri);
	static mpf_t cache_q;
	static int init = 0;

	if (!init)
	{
		init = 1;
		mpf_init (cache_q);
	}

	if (!mpf_eq(q,cache_q, prec*3.322))
	{
		fp_one_d_cache_clear (&re_diri);
		fp_one_d_cache_clear (&im_diri);
		mpf_set(cache_q,q);
	}

	if (prec <= fp_one_d_cache_check (&re_diri, k))
	{
		fp_one_d_cache_fetch (&re_diri, diri->re, k);
		fp_one_d_cache_fetch (&im_diri, diri->im, k);
		return;
	}
	
	diri_term_helper (diri, k, q, ess, prec);

	fp_one_d_cache_check (&im_diri, k);
	fp_one_d_cache_store (&re_diri, diri->re, k, prec);
	fp_one_d_cache_store (&im_diri, diri->im, k, prec);
}

/* ======================== END OF FILE ========================== */

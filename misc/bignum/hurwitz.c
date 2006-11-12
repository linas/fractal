/* 
 * hurwitz.c
 *
 * Compute the Hurwitz zeta function for arbitrary complex argument
 *
 * Linas Vepstas October 2006
 */

#include <math.h>
#include <stdio.h>

#include <gmp.h>
#include "mp-binomial.h"
#include "mp-cache.h"
#include "mp-complex.h"
#include "mp-misc.h"
#include "mp-trig.h"

/* =========================================================== */

/**
 * diri_term -- return (k+q)^{1-s}
 *
 * Values are cached, because they will be repeatedly called
 * for forward differencing.
 */
static void diri_term_helper (cpx_t *diri, int k, mpf_t q, cpx_t *ess, int prec)
{
	mpf_t kq, logkq, mag, pha;
	mpf_init (kq);
	mpf_init (logkq);
	mpf_init (mag);
	mpf_init (pha);

	mpf_add_ui (kq, q, k);
	
	fp_log (logkq, kq, prec);
	
	/* magnitude is exp((1-re(s)) *log(k+q)) */
	mpf_neg (mag, ess->re);
	mpf_add_ui (mag, mag, 1);
	mpf_mul (mag, mag, logkq);
	fp_exp (mag, mag, prec);

	/* phase is -im(s) *log(k+q)) */
	mpf_mul (pha, ess->im, logkq);
	mpf_neg (pha,pha);

	fp_cosine (diri->re, pha, prec);
	mpf_mul (diri->re, mag, diri->re);
	
	fp_sine (diri->im, pha, prec);
	mpf_mul (diri->im, mag, diri->im);
	
	mpf_clear(kq);
	mpf_clear(logkq);
	mpf_clear(mag);
	mpf_clear(pha);
}

static void diri_term (cpx_t *diri, int k, mpf_t q, cpx_t *ess, int prec)
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

/* =========================================================== */

/* Compute n'th forward difference of (k+q)^{1-s} */
static void forward_diff_diri (cpx_t *fin, int n, mpf_t q, cpx_t *ess, int prec)
{
	mpf_set_ui (fin->re, 0);
	mpf_set_ui (fin->im, 0);

	mpz_t ibin;
	mpz_init (ibin);
	
	mpf_t bin;
	mpf_init (bin);

	cpx_t diri;
	cpx_init (&diri);

	int k;
	for (k=0; k<=n; k++)
	{
		i_binomial (ibin, n,k);
		mpf_set_z (bin, ibin);

		diri_term (&diri, k, q, ess, prec);
		mpf_mul (diri.re, diri.re, bin);
		mpf_mul (diri.im, diri.im, bin);
		  
		if (0 == k%2)
		{
			cpx_add (fin, fin, &diri);
		}
		else
		{
			cpx_sub (fin, fin, &diri);
		}
	}

	mpz_clear (ibin);
	mpf_clear (bin);
	cpx_clear (&diri);
}

/* =========================================================== */
/* A brute-force summation using Hasse formula, 
 * for complex s, real q.
 *
 * Unfortunately, the convergence is slow on the critical strip.
 */

void hurwitz_zeta(cpx_t *zeta, cpx_t *ess, mpf_t q, int prec)
{
	int norder = prec;

	mpf_set_ui (zeta->re, 0);
	mpf_set_ui (zeta->im, 0);

	/* os contains value of 1/(s-1) */
	cpx_t os;
	cpx_init (&os);
	mpf_sub_ui (os.re,ess->re, 1);
	mpf_set (os.im, ess->im);
	cpx_recip (&os, &os);

	/* containes finite difference */
	cpx_t fd;
	cpx_init (&fd);

	mpf_t on;
	mpf_init (on);

	int n;
	for (n=0; n<norder; n++)
	{
		forward_diff_diri (&fd, n, q, ess, prec);
		
		mpf_set_ui (on, 1);
		mpf_div_ui (on, on, n+1);
		
		mpf_mul (fd.re, fd.re, on);
		mpf_mul (fd.im, fd.im, on);
printf ("duude %d ", n);
fp_prt (" ", fd.re);
		cpx_add (zeta, zeta, &fd);
cpx_mul (&fd, zeta, &os);
fp_prt ("   ", fd.re);
printf ("\n");
	}

	cpx_mul (zeta, zeta, &os);

	cpx_clear (&os);
	
	mpf_clear (on);
	cpx_clear (&fd);
} 

int main ()
{
	int prec = 580;

	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (3.3*prec+600);

	cpx_t ess, zeta;
	cpx_init (&ess);
	cpx_init (&zeta);

	mpf_set_d (ess.re, 0.5);
	mpf_set_d (ess.im, 4.0);
	
	mpf_t que;
	mpf_init (que);
	mpf_set_d (que,0.5);
	
	hurwitz_zeta (&zeta, &ess, que, prec);

	fp_prt ("its ", zeta.re);
	printf ("\n");

	return 0;
}

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
#include "mp-cache.h"
#include "mp-misc.h"
#include "mp-trig.h"

typedef struct {
	mpf_t re;
	mpf_t im;
} cpx_t;

static inline void cpx_init (cpx_t *z)
{
	mpf_init (z->re);
	mpf_init (z->im);
}

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
	
fp_prt ("its re ", diri->re);
printf ("\n");

fp_prt ("its im ", diri->im);
printf ("\n");

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

#if 0
/* A brute-force summation using Hasse formula, 
 * for complex s, real q.
 *
 * Unfortunately, the convergence is slow on the critical strip,
 * and double precision is not enough to do anything useful here.
 */

void hurwitz_zeta (long double *phre, long double *phim, double sre, double sim, long double q)
{
	int norder = 80;

	long double hre = 0.0;
	long double him = 0.0;
	int n;
	for (n=0; n<norder; n++)
	{
		long double rs=0.0L;
		long double is=0.0L;
		long double cyn = 1.0L;
		for (k=0; k<=n; k++)
		{
			long double bin = cyn*binomial (n,k);
			rs += bin * refd[k];
			is += bin * imfd[k];
			cyn = -cyn;
		}

		long double on = 1.0L/(n+1.0L);
		hre += on * rs;
		him += on * is;
		printf ("its %d \t%Lg \thre=%Lg \t%Lg \thim=%Lg\n", n, rs, hre,is, him);
	}
} 
#endif

int main ()
{

#if 0
	long double q;
	long double sre, sim;
	long double hre, him;

	sre = 0.5;
	sim = 2.0;
	q = 0.3;
	hurwitz_zeta (&hre, &him, sre,sim, q);
#endif

	cpx_t ess, diri;
	cpx_init (&ess);
	cpx_init (&diri);

	mpf_set_ui (ess.re, 3);
	
	mpf_t que;
	mpf_init (que);
	mpf_set_d (que,0.5);
	
	diri_term (&diri, 2, que, &ess, 50);

	fp_prt ("its ", diri.re);
	printf ("\n");

	return 0;
}

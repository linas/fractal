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

typedef struct {
	mpf_t re;
	mpf_t im;
} cpx_t;

static inline void cpx_init (cpx_t z)
{
	mpf_init (z.re);
	mpf_init (z.im);
}

/**
 * diri_term -- return (k+q)^{-s}
 *
 * Values are cached, because they will be repeatedly called
 * for forward differencing.
 */
static void diri_term (cpx_t diri, int k, double q)
{
	DECLARE_FP_CACHE (re_diri);
	DECLARE_FP_CACHE (im_diri);

	int k;
	for (k=0; k<norder; k++)
	{
		long double logkq = logl(k+q);
		long double mag = expl((1.0L-sre) * logkq);
		refd[k] = mag * cosl (sim*logkq);
		imfd[k] = mag * sinl (sim*logkq);

		// printf ("its %d \t%Lg \t%Lg\n", k, refd[k], imfd[k]);
	}
}



/* A brute-force summation using Hasse formula, 
 * for complex s, real q.
 *
 * Unfortunately, the convergence is slow on the critical strip,
 * and double precision is not enough to do anything useful here.
 */

void hurwitz_zeta (long double *phre, long double *phim, double sre, double sim, long double q)
{
	int norder = 80;

	/* arrays storing values to be forward-difference'd */
	DECLARE_FP_CACHE (refd);
	DECLARE_FP_CACHE (imfd);

	int k;
	for (k=0; k<norder; k++)
	{
		long double logkq = logl(k+q);
		long double mag = expl((1.0L-sre) * logkq);
		refd[k] = mag * cosl (sim*logkq);
		imfd[k] = mag * sinl (sim*logkq);

		// printf ("its %d \t%Lg \t%Lg\n", k, refd[k], imfd[k]);
	}

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

main ()
{
	long double q;
	long double sre, sim;
	long double hre, him;

	sre = 0.5;
	sim = 2.0;
	q = 0.3;
	hurwitz_zeta (&hre, &him, sre,sim, q);
}

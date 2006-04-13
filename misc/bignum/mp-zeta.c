/*
 * mp_zeta.c
 *
 * High-precison Riemann zeta function, using the 
 * Gnu Multiple-precision library.
 *
 * Also, high-precision values of the series a_n 
 * 
 * Linas Vepstas July 2005
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mp_zeta.h"

void i_prt (char * str, mpz_t val)
{
	printf (str);
	mpz_out_str (stdout, 10, val);
	printf ("\n");
}

void fp_prt (char * str, mpf_t val)
{
	printf (str);
	mpf_out_str (stdout, 10, 60, val);
	printf ("\n");
}

/* ======================================================================= */
/* Cache management */

typedef struct {
	unsigned int nmax;
	mpz_t *cache;
	char *ticky;
	short disabled;
} i_cache;


#define DECLARE_I_CACHE(name)         \
	static i_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, .disabled = 0}

/** i_one_d_cache_check() -- check if mpz_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int i_one_d_cache_check (i_cache *c, unsigned int n)
{
	if (c->disabled) return 0;
	if (n > c->nmax)
	{
		unsigned int newsize = n+1;
		c->cache = (mpz_t *) realloc (c->cache, newsize * sizeof (mpz_t));
		c->ticky = (char *) realloc (c->ticky, newsize * sizeof (char));

		unsigned int en;
		for (en=c->nmax; en <newsize; en++)
		{
			mpz_init (c->cache[en]);
			c->ticky[en] = 0;
		}
		c->nmax = n;
		return 0;
	}

	return (c->ticky[n]);
}

/** 
 * i_one_d_cache_fetch - fetch value from cache
 */
void i_one_d_cache_fetch (i_cache *c, mpz_t val, unsigned int n)
{
	if (c->disabled) return;
	mpz_set (val, c->cache[n]);
}

/**
 * i_one_d_cache_store - store value in cache
 */
void i_one_d_cache_store (i_cache *c, mpz_t val, unsigned int n)
{
	if (c->disabled) return;
	mpz_set (c->cache[n], val);
	c->ticky[n] = 1;
}

/* ======================================================================= */
/** i_triangle_cache_check() -- check if mpz_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a trianglular cache layout (two indecies)
 *  with 0 <= k <=n
 */
int i_triangle_cache_check (i_cache *c, unsigned int n, unsigned int k)
{
	if (n> c->nmax)
	{
		unsigned int newsize = (n+1)*(n+2)/2;
		c->cache = (mpz_t *) realloc (c->cache, newsize * sizeof (mpz_t));
		c->ticky = (char *) realloc (c->ticky, newsize * sizeof (char));

		unsigned int en;
		for (en=c->nmax+1; en <=n; en++)
		{
			unsigned int j;
			unsigned int idx = en * (en+1) /2;
			for (j=0; j<=en; j++)
			{
				mpz_init (c->cache[idx+j]);
				c->ticky[idx+j]=0;
			}
		}
		c->nmax = n;
		return 0;
	}
	unsigned int idx = n * (n+1) /2 ;
	return c->ticky[idx+k];
}

/** 
 * i_triangle_cache_fetch - fetch value from cache
 */
void i_triangle_cache_fetch (i_cache *c, mpz_t val, unsigned int n, unsigned int k)
{
	unsigned int idx = n * (n+1) /2 ;
	mpz_set (val, c->cache[idx+k]);
}

/**
 * i_triangle_cache_store - store value in cache
 */
void i_triangle_cache_store (i_cache *c, mpz_t val, unsigned int n, unsigned int k)
{
	unsigned int idx = n * (n+1) /2 ;
	mpz_set (c->cache[idx+k], val);
	c->ticky[idx+k] = 1;
}

/* ======================================================================= */
/* Cache management */
/* pure cut-n-paste of he integer variant */

typedef struct {
	unsigned int nmax;
	mpq_t *cache;
	char *ticky;
} q_cache;

#define DECLARE_Q_CACHE(name)         \
	static q_cache name = {.nmax=0, .cache=NULL, .ticky=NULL}

/** q_one_d_cache_check() -- check if mpq_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int q_one_d_cache_check (q_cache *c, unsigned int n)
{
	if (n > c->nmax)
	{
		unsigned int newsize = n+1;
		c->cache = (mpq_t *) realloc (c->cache, newsize * sizeof (mpq_t));
		c->ticky = (char *) realloc (c->ticky, newsize * sizeof (char));

		unsigned int en;
		for (en=c->nmax; en <newsize; en++)
		{
			mpq_init (c->cache[en]);
			c->ticky[en] = 0;
		}
		c->nmax = n;
		return 0;
	}

	return (c->ticky[n]);
}

/** 
 * q_one_d_cache_fetch - fetch value from cache
 */
void q_one_d_cache_fetch (q_cache *c, mpq_t val, unsigned int n)
{
	mpq_set (val, c->cache[n]);
}

/**
 * q_one_d_cache_store - store value in cache
 */
void q_one_d_cache_store (q_cache *c, mpq_t val, unsigned int n)
{
	mpq_set (c->cache[n], val);
	c->ticky[n] = 1;
}

/* ======================================================================= */
/* Cache management */
/* Almost a cut-n-paste of above, but using fp instead */

typedef struct {
	unsigned int nmax;
	mpf_t *cache;
	int *precision; /* base-10 precision of cached value */
} fp_cache;


#define DECLARE_FP_CACHE(name)         \
	static fp_cache name = {.nmax=0, .cache=NULL, .precision=NULL}

/** fp_one_d_cache_check() -- check if mpf_t value is in the cache
 *  If there is a cached value, this returns the precision of the 
 *  value in the cache; else it returns zero.
 *  This assumes a 1-dimensional cache layout (simple array)
 */
int fp_one_d_cache_check (fp_cache *c, unsigned int n)
{
	if (n > c->nmax)
	{
		unsigned int newsize = n+1;
		c->cache = (mpf_t *) realloc (c->cache, newsize * sizeof (mpf_t));
		c->precision = (int *) realloc (c->precision, newsize * sizeof (int));

		unsigned int en;
		for (en=c->nmax; en <newsize; en++)
		{
			mpf_init (c->cache[en]);
			c->precision[en] = 0;
		}
		c->nmax = n;
		return 0;
	}

	return (c->precision[n]);
}

/** 
 * fp_one_d_cache_fetch - fetch value from cache
 */
void fp_one_d_cache_fetch (fp_cache *c, mpf_t val, unsigned int n)
{
	mpf_set (val, c->cache[n]);
}

/**
 * fp_one_d_cache_store - store value in cache
 */
void fp_one_d_cache_store (fp_cache *c, mpf_t val, unsigned int n, int prec)
{
	mpf_set (c->cache[n], val);
	c->precision[n] = prec;
}

void fp_one_d_cache_clear (fp_cache *c)
{
	int i;
	for (i=0; i<c->nmax; i++)
	{
		c->precision[i] = 0;
	}
}

/* ======================================================================= */
/** fp_triangle_cache_check() -- check if mpf_t value is in the cache
 *  If there is a cached value, this returns the precision of the 
 *  value in the cache; else it returns zero.
 *  This assumes a trianglular cache layout (two indecies)
 *  with 0 <= k <=n
 */
int fp_triangle_cache_check (fp_cache *c, unsigned int n, unsigned int k)
{
	if (n> c->nmax)
	{
		unsigned int newsize = (n+1)*(n+2)/2;
		c->cache = (mpf_t *) realloc (c->cache, newsize * sizeof (mpf_t));
		c->precision = (int *) realloc (c->precision, newsize * sizeof (int));

		unsigned int en;
		for (en=c->nmax+1; en <=n; en++)
		{
			unsigned int j;
			unsigned int idx = en * (en+1) /2 ;
			for (j=0; j<=en; j++)
			{
				mpf_init (c->cache[idx+j]);
				c->precision[idx+j]=0;
			}
		}
		c->nmax = n;
		return 0;
	}
	unsigned int idx = n * (n+1) /2 ;
	return c->precision[idx+k];
}

/** 
 * fp_triangle_cache_fetch - fetch value from cache
 */
void fp_triangle_cache_fetch (fp_cache *c, mpf_t val, unsigned int n, unsigned int k)
{
	unsigned int idx = n * (n+1) /2 ;
	mpf_set (val, c->cache[idx+k]);
}

/**
 * fp_triangle_cache_store - store value in cache
 */
void fp_triangle_cache_store (fp_cache *c, mpf_t val, 
					 unsigned int n, unsigned int k, int prec)
{
	unsigned int idx = n * (n+1) /2 ;
	mpf_set (c->cache[idx+k], val);
	c->precision[idx+k] = prec;
}

/* ======================================================================= */
/* i_poch_rising
 * rising pochhammer symbol, for integer values.
 *
 * Brute force, simple.
 */

void i_poch_rising (mpz_t poch, unsigned int k, unsigned int n)
{
	mpz_t acc;
	mpz_init (acc);

	mpz_set_ui (poch, 1);
	unsigned int i;
	for (i=0; i<n; i++)
	{
		mpz_mul_ui (acc, poch, i+k);
		mpz_set (poch, acc);
	}

	mpz_clear (acc);
}

/** 
 * i_factorial -- the factorial
 */
#ifdef USE_LOCAL_FACTORIAL
void i_factorial (mpz_t fact, unsigned int n)
{
	DECLARE_I_CACHE (cache);

	if (1 >= n)
	{
		mpz_set_ui (fact, 1);
		return;
	}
	int hit = i_one_d_cache_check (&cache, n);
	if (hit)
	{
		i_one_d_cache_fetch (&cache, fact, n);
	}
	else
	{
		i_poch_rising (fact, 1, n);
		i_one_d_cache_store (&cache, fact, n);
	}
}
#else
#define i_factorial mpz_fac_ui
#endif /* USE_LOCAL_FACTORIAL */

/* ====================================================================== */
/* i_binomial
 * Binomial coefficient (n k)
 */

#ifdef USE_LOCAL_BINOMIAL
void i_binomial_compute (mpz_t bin, unsigned int n, unsigned int k)
{
	mpz_t top, bot;

	if (2*k < n) k = n-k;

	mpz_init (top);
	mpz_init (bot);
	i_poch_rising (top, k+1, n-k);
	i_factorial (bot, n-k); 

	mpz_divexact (bin, top, bot);
	
	mpz_clear (top);
	mpz_clear (bot);
}

/**
 * i_binomial - return the binomial coefficient
 * Uses a cached value if avalable.
 */ 
void i_binomial (mpz_t bin, unsigned int n, unsigned int k)
{
	DECLARE_I_CACHE (cache);

	if (1 >= n)
	{
		mpz_set_ui (bin, 1);
		return;
	}

	if (2*k < n) k = n-k;
	int hit = i_triangle_cache_check (&cache, n, k);
	if (hit)
	{
		i_triangle_cache_fetch (&cache, bin, n, k);
	}
	else
	{
		i_binomial_compute (bin, n, k);
		i_triangle_cache_store (&cache, bin, n, k);
	}
}
#else 
#define i_binomial mpz_bin_uiui

#endif /* USE_LOCAL_BINOMIAL */

/* ======================================================================= */
/* stirling_first - Stirling Numbers of the First kind, 
 * normalized so that they are all positive.
 * Uses dynamically-sized cache.
 */
void i_stirling_first (mpz_t s, unsigned int n, unsigned int k)
{
	DECLARE_I_CACHE (cache);

	/* Trivial case (not in the cache) */
	if (0==k)
	{
		if (0==n) 
		{ 
			mpz_set_ui (s, 1);
		}
		else
		{
			mpz_set_ui (s, 0);
		}
		return;
	}

	if (n<k)
	{
		mpz_set_ui (s, 0);
		return;
	}

	if (n==k)
	{
		mpz_set_ui (s, 1);
		return;
	}

	/* Pull value from cache if it is there */
	int hit = i_triangle_cache_check (&cache, n, k);
	if (hit)
	{
		i_triangle_cache_fetch (&cache, s, n, k);
		return;
	}
	
	/* Use recursion to get new value */
	/* s = stir(n-1, k-1) + (n-1) * stir(n-1, k) */
	unsigned int i;
	mpz_t skm, sk, en;
	mpz_init (skm);
	mpz_init (sk);
	mpz_init (en);
	mpz_set_ui (skm, 0);
	mpz_set_ui (en, n-1);
	for (i=1; i<=n; i++)
	{
		i_stirling_first (sk, n-1, i);
		mpz_mul (s, en, sk);
		mpz_add (s, s, skm);
		i_triangle_cache_store (&cache, s, n, i);
		mpz_set (skm, sk);
	}
	mpz_clear (skm);
	mpz_clear (sk);
	mpz_clear (en);

	i_triangle_cache_fetch (&cache, s, n, k);
}

/* ======================================================================= */
/* A funny off-by-one sum of stirling and binomial */

void i_stirbin_sum_compute (mpz_t s, unsigned int n, unsigned int m)
{
	unsigned int k;
	
	mpz_t term, stir, bin;
	mpz_init (term);
	mpz_init (stir);
	mpz_init (bin);
	mpz_set_ui (s, 0);
	for (k=m; k<=n; k++)
	{
		i_stirling_first (stir, n, k);
		i_binomial (bin, k,m);
		mpz_mul (term, bin, stir);
		if (k%2)
		{
			mpz_sub (s, s, term);
		}
		else
		{
			mpz_add (s, s, term);
		}
	}
	mpz_clear (term);
	mpz_clear (stir);
	mpz_clear (bin);
}

void i_stirbin_sum (mpz_t s, unsigned int n, unsigned int m)
{
	DECLARE_I_CACHE (cache);

	if (0 >= n)
	{
		mpz_set_ui (s, 1);
		return;
	}

	int hit = i_triangle_cache_check (&cache, n, m);
	if (hit)
	{
		i_triangle_cache_fetch (&cache, s, n, m);
	}
	else
	{
		i_stirbin_sum_compute (s, n, m);
		i_triangle_cache_store (&cache, s, n, m);
	}
}

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
/* binomial transform of power sum */

void fp_bin_xform_pow_compute (mpf_t bxp, unsigned int n, unsigned int s)
{
	mpz_t bin;
	mpz_init (bin);

	mpf_t vp, term;
	mpf_init (vp);
	mpf_init (term);
	
	mpf_set_ui (bxp, 0);
	int k;
	for (k=0; k<=n; k++)
	{
		i_binomial (bin, n, k);
		mpf_set_z (term, bin);
		fp_inv_pow (vp, k+1, s);
		mpf_mul (term, term, vp);

		if (k%2)
		{
			mpf_sub (bxp, bxp, term);
		}
		else
		{
			mpf_add (bxp, bxp, term);
		}
	}
	mpz_clear (bin);
	mpf_clear (vp);
	mpf_clear (term);
}

void fp_bin_xform_pow (mpf_t bxp, unsigned int n, unsigned int s)
{
	DECLARE_FP_CACHE (cache);
	if (0 == n)
	{
		mpf_set_ui (bxp, 1); 
		return;
	}
	int hit = fp_triangle_cache_check (&cache, n+s, s);
	if (hit)
	{
		fp_triangle_cache_fetch (&cache, bxp, n+s, s);
	}
	else
	{
		fp_bin_xform_pow_compute (bxp, n, s);
		fp_triangle_cache_store (&cache, bxp, n+s, s, 1);
	}
}

/* ======================================================================= */
/* Helmut Hasse zeta -- use the globally convergent Helmut Hasse zeta
 * series for computing the value of the zeta ... 
 */

void fp_hasse_zeta_compute (mpf_t zeta, unsigned int s, int prec)
{
	// This gets the decimal pecision just right!
	// This works because the bin_xform_pow is always of order 1.
	int nmax = (3.321928*prec+3);
	int n;
	
	mpf_t twon, term;
	mpf_init (twon);
	mpf_init (term);
	
	mpf_set_ui (twon, 1);
	mpf_div_ui (twon, twon, 2);
	
	mpf_set_ui (zeta, 0);
	for (n=0; n<nmax; n++)
	{
		fp_bin_xform_pow (term, n, s);
		mpf_mul (term, term, twon);
		mpf_add (zeta, zeta, term);

		mpf_div_ui (twon, twon, 2);
	}
	mpf_set_ui (twon, 1);
	mpf_div_2exp (twon, twon, s-1);
	
	mpf_set_ui (term, 1);
	mpf_sub (term, term, twon);
	
	mpf_div (zeta, zeta, term);

	mpf_clear (twon);
	mpf_clear (term);
}

void fp_hasse_zeta (mpf_t zeta, unsigned int s, int prec)
{
	DECLARE_FP_CACHE (cache);
	if (1 >= s)
	{
		mpf_set_ui (zeta, 0); 
		return;
	}
	int have_prec = fp_one_d_cache_check (&cache, s);
	if (have_prec >= prec)
	{
		fp_one_d_cache_fetch (&cache, zeta, s);
	}
	else
	{
		fp_hasse_zeta_compute (zeta, s, prec);
		fp_one_d_cache_store (&cache, zeta, s, prec);
	}
}

/* ======================================================================= */
/** 
 * fp_harmonic -- The harmonic number
 */
void fp_harmonic (mpf_t harm, unsigned int n)
{
	DECLARE_FP_CACHE (cache);

	if (1 >= n)
	{
		mpf_set_ui (harm, 1);
		return;
	}
	int hit = fp_one_d_cache_check (&cache, n);
	if (hit)
	{
		fp_one_d_cache_fetch (&cache, harm, n);
		return;
	}
	
	int istart = n-1;
	hit = fp_one_d_cache_check (&cache, istart);
	while (0 == hit && istart>1)
	{
		istart--;
		hit = fp_one_d_cache_check (&cache, istart);
	}

	int i;
	fp_harmonic (harm, istart);

	mpf_t term;
	mpf_init (term);
	for (i=istart+1; i<=n; i++)
	{
		mpf_set_ui (term, 1);
		mpf_div_ui (term, term, i);
		mpf_add (harm, harm, term);
		fp_one_d_cache_store (&cache, harm, i, 1);
	}
	mpf_clear (term);
}

/* ======================================================================= */
/* fp_poch_rising
 * rising pochhammer symbol (x)_n, for real values of x and integer n.
 *
 * Brute force, simple.
 */

void fp_poch_rising (mpf_t poch, double x, unsigned int n)
{
	mpf_t acc, term;
	
	mpf_init (acc);
	mpf_init (term);

	mpf_set_ui (poch, 1);
	unsigned int i;
	for (i=0; i<n; i++)
	{
		mpf_set_d (term, x+i);
		mpf_mul (acc, poch, term);
		mpf_set (poch, acc);
	}

	mpf_clear (acc);
	mpf_clear (term);
}

/* ======================================================================= */
/* c_poch_rising
 * rising pochhammer symbol (s)_n, for complex s and integer n.
 *
 * Brute force, simple.
 */

void c_poch_rising (mpf_t re_poch, mpf_t im_poch, double re_s, double im_s, unsigned int n)
{
	mpf_t racc, iacc, atmp, btmp, re_term, im_term;
	
	mpf_init (racc);
	mpf_init (iacc);
	mpf_init (atmp);
	mpf_init (btmp);
	mpf_init (re_term);
	mpf_init (im_term);

	mpf_set_ui (re_poch, 1);
	mpf_set_ui (im_poch, 0);
	unsigned int i;
	for (i=0; i<n; i++)
	{
		mpf_set_d (re_term, re_s+i);
		mpf_set_d (im_term, im_s);

		mpf_mul (atmp, re_poch, re_term);
		mpf_mul (btmp, im_poch, im_term);
		mpf_sub (racc, atmp, btmp);

		mpf_mul (atmp, re_poch, im_term);
		mpf_mul (btmp, im_poch, re_term);
		mpf_add (iacc, atmp, btmp);

		mpf_set (re_poch, racc);
		mpf_set (im_poch, iacc);
	}

	mpf_clear (racc);
	mpf_clear (iacc);
	mpf_clear (atmp);
	mpf_clear (btmp);
	mpf_clear (re_term);
	mpf_clear (im_term);
}

/* ======================================================================= */
/* fp_binomial
 * Binomial coefficient 
 */

void fp_binomial (mpf_t bin, double s, unsigned int k)
{
	mpf_t top, bot;
	mpz_t fac;

	mpf_init (top);
	mpf_init (bot);
	mpz_init (fac);
	fp_poch_rising (top, s-k+1, k);
	i_factorial (fac, k); 
	mpf_set_z (bot, fac);

	mpf_div (bin, top, bot);
	
	mpf_clear (top);
	mpf_clear (bot);
	mpz_clear (fac);
}

/* ======================================================================== */
/* c_binomial
 * Complex binomial coefficient
 */

void c_binomial (mpf_t re_bin, mpf_t im_bin, double re_s, double im_s, unsigned int k)
{
	mpf_t retop, imtop, bot;
	mpz_t fac;

	mpf_init (retop);
	mpf_init (imtop);
	mpf_init (bot);
	mpz_init (fac);
	c_poch_rising (retop, imtop, re_s-k+1, im_s, k);
	i_factorial (fac, k); 
	mpf_set_z (bot, fac);

	mpf_div (re_bin, retop, bot);
	mpf_div (im_bin, imtop, bot);
	
	mpf_clear (retop);
	mpf_clear (imtop);
	mpf_clear (bot);
	mpz_clear (fac);
}

/* ======================================================================= */
/* Bernoulli number as a rational */

void q_bernoulli (mpq_t bern, int n)
{
	DECLARE_Q_CACHE (cache);

	if (0>n) return;
	if (0==n) {	mpq_set_ui (bern, 1,1); return; }
	if (1==n) { mpq_set_si (bern, -1,2); return; }

	/* All other odd n's are zero */
	if (n%2) { mpq_set_ui (bern, 0, 1);  return; }

	int hn = n/2;

	int hit = q_one_d_cache_check (&cache, hn);
	if (hit)
	{
		q_one_d_cache_fetch (&cache, bern, hn);
		return;
	}
	
	/* Not found in cache, will have to compute */
	mpz_t binom;
	mpz_init (binom);

	mpq_t term, tmp;
	mpq_init (term);
	mpq_init (tmp);

	mpq_set_si (bern, 1-n, 2);
	
	int i;
	for (i=1; i<hn; i++)
	{
		int k = 2*i;
		i_binomial (binom, n+1, k);

		mpq_set_z (tmp, binom);
		q_bernoulli (term, k);
		mpq_mul (term, term, tmp);
		mpq_add (bern, bern, term);
	}

	mpq_set_si (tmp, -1, n+1);
	mpq_mul (bern, bern, tmp);

	mpz_clear (binom);
	mpq_clear (term);
	mpq_clear (tmp);

	q_one_d_cache_store (&cache, bern, hn);
}

/* ======================================================================= */
/**
 * fp_exp -  Floating point exponential
 * Implemented using a brute-force, very simple algo, with 
 * no attempts at optimization. Also, does not assume any 
 * precomputed constants.
 */

void fp_exp (mpf_t ex, mpf_t z, unsigned int prec)
{
	mpf_t z_n, fact, term;

	mpf_init (z_n);
	mpf_init (fact);
	mpf_init (term);

	mpf_set_ui (ex, 1);
	mpf_set_ui (fact, 1);
	mpf_set (z_n, z);
	
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
		
		/* don't go no father than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n++;
		mpf_mul (z_n, z_n, z);
		mpf_mul_ui (fact, fact, n);
	}
	
	mpf_clear (z_n);
	mpf_clear (fact);
	mpf_clear (term);

	mpf_clear (one);
	mpf_clear (maxterm);
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
	mpf_t z_n, zsq, term;

	mpf_init (z_n);
	mpf_init (zsq);
	mpf_init (term);

	mpf_set_ui (atn, 0);
	mpf_mul (zsq, z, z);
	mpf_set (z_n, z);
	
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
		mpf_div_ui (term, z_n, n);
		if (n%2)
		{
			mpf_add (atn, atn, term);
		}
		else
		{
			mpf_sub (atn, atn, term);
		}
		
		/* don't go no father than this */
		mpf_abs (term, term);
		if (mpf_cmp (term, maxterm) < 0) break;
		
		n += 2;
		mpf_mul (z_n, z_n, zsq);
	}
	
	mpf_clear (z_n);
	mpf_clear (zsq);
	mpf_clear (term);

	mpf_clear (one);
	mpf_clear (maxterm);
}

/* ======================================================================= */
/* fp_euler
 * return Euler-Mascheroni const
 */
void fp_euler_mascheroni (mpf_t gam)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		// char * g="0.5772156649015328606065120900824024310421593359399235988057672348848677267776646709369470632917467495e0";
	char * g="0.57721566490153286060651209008240243104215933593992359880576723488486772677766 467093694706329174674951463144724980708248096050401448654283622417399764492353 625350033374293733773767394279259525824709491600873520394816567085323315177661 152862119950150798479374508570574002992135478614669402960432542151905877553526 733139925401296742051375413954911168510280798423487758720503843109399736137255 306088933126760017247953783675927135157722610273492913940798430103417771778088 154957066107501016191663340152278935867965497252036212879226555953669628176388";
	
		mpf_set_str (e, g, 10);
	}
	mpf_set (gam, e);
}

void fp_pi_string (mpf_t pi)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		char *p = "3.1415926535 8979323846 2643383279 5028841971 6939937510"
			   "5820974944 5923078164 0628620899 8628034825 3421170679"
				 "8214808651 3282306647 0938446095 5058223172 5359408128"
				  "4811174502 8410270193 8521105559 6446229489 5493038196"
				   "4428810975 6659334461 2847564823 3786783165 2712019091"
					 "4564856692 3460348610 4543266482 1339360726 0249141273"
					  "7245870066 0631558817 4881520920 9628292540 9171536436"
					   "7892590360 0113305305 4882046652 1384146951 9415116094"
						 "3305727036 5759591953 0921861173 8193261179 3105118548"
						  "0744623799 6274956735 1885752724 8912279381 8301194912"
						   "9833673362 4406566430 8602139494 6395224737 1907021798"
							 "6094370277 0539217176 2931767523 8467481846 7669405132"
							  "0005681271 4526356082 7785771342 7577896091 7363717872"
							   "1468440901 2249534301 4654958537 1050792279 6892589235"
								 "4201995611 2129021960 8640344181 5981362977 4771309960"
								  "5187072113 4999999837 2978049951 0597317328 1609631859"
								   "5024459455 3469083026 4252230825 3344685035 2619311881";
		mpf_set_str (e, p, 10);
	}
	mpf_set (pi, e);
}

void fp_pi (mpf_t pi, int prec)
{
	static int precision=0;
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

/* return e^pi */
void fp_e_pi (mpf_t e_pi)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		char *p= "23.140692632779269005729086367948547380266106242600211993445046409524342350690 452783516971997067549219675952704801087773144428044414693835844717445879609849 365327965863669242230268991013741764684401410395183868477243068059588162449844 491430966778413671631963414784038216511287637731470347353833162821294047891936 224820221006032065443362736557271823744989618858059591684872645479013397834026 595101499643792422968160799565381423536206957600770590460899883002254304871211 791300849327379580729427301931042601691939325853203428968661895283290521711157 185185506802254197204566370865568386830544799278170407497768540367556534957218 867882563994384718224585889428535247260568210271076018491534518468064887386774 439630514005169440540665265430968869063937315359837311042174433023967896690035";
		mpf_set_str (e, p, 10);
	}
	mpf_set (e_pi, e);
}

void fp_zeta2 (mpf_t zeta)
{
	static int inited=0;
	static mpf_t z;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (z);

		mpf_t pi, pisq;
		mpf_init (pi);
		mpf_init (pisq);
		
		fp_pi (pi);
		mpf_mul (pisq, pi, pi);
		mpf_div_ui (z, pisq, 6);
	
		mpf_clear (pi);
		mpf_clear (pisq);
	}
	mpf_set (zeta, z);
}

void fp_zeta3 (mpf_t zeta)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
		// char * g="1.202056903159594285399738161511449990764986292";
		char * g="1.2020569031595942853997381615114499907649862923404988817922715553418382057863 130901864558736093352581461991577952607194184919959986732832137763968372079001 614539417829493600667191915755222424942439615639096641032911590957809655146512 799184051057152559880154371097811020398275325667876035223369849416618110570147 157786394997375237852779370309560257018531827900030765471075630488433208697115";
		mpf_set_str (e, g, 10);
	}
	mpf_set (zeta, e);
}

void fp_zeta5 (mpf_t zeta)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
		char * g="1.036927755143369926331365486457034168057080919501912811974192677 9038035897862814845600431065571333363796203414665566090428009617 7915597084183511072180087644866286337180353598363962365128888981 3352767752398275032022436845766444665958115993917977745039244643 9196666159664016205325205021519226713512567859748692860197447984 3200672681297530919900774656558601526573730037561532683149897971 9350398378581319922884886425335104251602510849904346402941172432 7576341508162332245618649927144272264614113007580868316916497918";

		mpf_set_str (e, g, 10);
	}
	mpf_set (zeta, e);
}

void fp_zeta7 (mpf_t zeta)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
		char * g="1.008349277381922826839797549849796759599863560565238706417283136 5716014783173557353460969689138513239689614536514910748872867774 1984033544031579830103398456212106946358524390658335396467699756 7696691427804314333947495215378902800259045551979353108370084210 7329399046107085641235605890622599776098694754076320000481632951 2586769250630734413632555601360305007373302413187037951026624779 3954650225467042015510405582224239250510868837727077426002177100 0195455778989836046745406121952650765461161356548679150080858554";
		mpf_set_str (e, g, 10);
	}
	mpf_set (zeta, e);
}

void fp_zeta9 (mpf_t zeta)
{
	static int inited=0;
	static mpf_t e;

	if (0 == inited)
	{
		inited = 1;
		mpf_init (e);

		// http://www.worldwideschool.org/library/books/sci/math/MiscellaneousMathematicalConstants/chap97.html
		char * g="1.002008392826082214417852769232412060485605851394888756548596615 9097850533902583989503930691271695861574086047658470602614253739 7072243015306913249876425109092948687676545396979415407826022964 1544836250668629056707364521601531424421326337598815558052591454 0848901539527747456133451028740613274660692763390016294270864220 1123162209241265753326205462293215454665179945038662778223564776 1660330281492364570399301119383985017167926002064923069795850945 8457966548540026945118759481561430375776154443343398399851419383";

		mpf_set_str (e, g, 10);
	}
	mpf_set (zeta, e);
}

/* ======================================================================= */

void fp_zeta_even_ui (mpf_t zeta, unsigned int n, unsigned int div)
{
	mpf_t pi, pip;
	mpf_init (pi);
	mpf_init (pip);
	
	fp_pi (pi);
	mpf_pow_ui (pip, pi, n);
	mpf_div_ui (zeta, pip, div);

	mpf_clear (pi);
	mpf_clear (pip);
}

void fp_zeta_even_str (mpf_t zeta, unsigned int n, char * snum, char * sdenom)
{
	mpf_t pi, pip, num, denom;
	mpf_init (pi);
	mpf_init (pip);
	mpf_init (num);
	mpf_init (denom);

	mpf_set_str (num,snum, 10);
	mpf_set_str (denom, sdenom, 10);
	
	fp_pi (pi);
	mpf_pow_ui (pip, pi, n);
	mpf_mul(pi, pip, num);
	mpf_div(zeta, pi, denom);

	mpf_clear (pi);
	mpf_clear (pip);
	mpf_clear (num);
	mpf_clear (denom);
}

/* ======================================================================= */
/* Compute and return the "exact" result for the zeta function for 
 * any value of even n. This is obtained by recursievly computing
 * the Bernoulli numbers, and multiplying by an appropriate factor
 * of pi and factorial. 
 */
void fp_zeta_even (mpf_t zeta, unsigned int n)
{
	mpq_t bern, b2, bb;
	mpq_init (bern);
	mpq_init (b2);
	mpq_init (bb);

	q_bernoulli (bern, n);
	mpq_set_ui (bb, 1, 2);
	mpq_mul (b2, bern, bb);

	/* divide by factorial */
	mpz_t fact;
	mpz_init (fact);
	i_factorial (fact, n);
	mpq_set_z (bern, fact);
	mpq_div (bb, b2, bern);

	/* fix the sign */
	if (0==n%4) mpq_neg (bb, bb);
	
	mpf_t pi, pip;
	mpf_init (pi);
	mpf_init (pip);
	
	fp_pi (pi);
	mpf_mul_ui (pi, pi, 2);
	mpf_pow_ui (pip, pi, n);

	mpf_set_q (pi, bb);
	mpf_mul (zeta,pi, pip);

	mpf_clear (pi);
	mpf_clear (pip);

	mpq_clear (bern);
	mpq_clear (b2);
	mpq_clear (bb);

	mpz_clear (fact);
}

/* ======================================================================= */
/* Return sum_n (n^k (e^{\pi k} \pm 1)^{-1}
 * The Simon Plouffe Ramanujan inspired thingy
 */
static void fp_ess (mpf_t ess_plus, mpf_t ess_minus, unsigned int k, unsigned int prec)
{
	mpf_t e_pi, en, enp, epip, eppos, epneg, term, oterm, acc;

	mpf_init (e_pi);
	mpf_init (en);
	mpf_init (enp);
	mpf_init (epip);
	mpf_init (eppos);
	mpf_init (epneg);
	mpf_init (term);
	mpf_init (oterm);
	mpf_init (acc);

	fp_e_pi (e_pi);
	mpf_set_ui (ess_plus, 0);
	mpf_set_ui (ess_minus, 0);

	// double mex = ((double) prec) * log (10.0) / log(2.0);
	double mex = ((double) prec) * 3.321928095;
	unsigned int imax = (unsigned int) (mex +1.0);
	mpf_t maxterm, one;
	mpf_init (maxterm);
	mpf_init (one);
	mpf_set_ui (one, 1);
	mpf_mul_2exp (maxterm, one, imax);
	
	int n;
	for (n=1; n<1000000000; n++)
	{
		mpf_set_ui (en, n);
		mpf_pow_ui (enp, en, k);
		mpf_pow_ui (epip, e_pi, 2*n);
		
		mpf_add_ui (eppos, epip, 1);
		mpf_sub_ui (epneg, epip, 1);

		mpf_mul (term, enp, eppos);
		mpf_ui_div (oterm, 1, term);
		mpf_add (acc, ess_plus, oterm);
		mpf_set (ess_plus, acc);

		mpf_mul (term, enp, epneg);
		mpf_ui_div (oterm, 1, term);
		mpf_add (acc, ess_minus, oterm);
		mpf_set (ess_minus, acc);

		/* don't go no father than this */
		if (mpf_cmp (term, maxterm) > 0) break;
	}

	mpf_clear (e_pi);
	mpf_clear (en);
	mpf_clear (enp);
	mpf_clear (epip);
	mpf_clear (eppos);
	mpf_clear (epneg);
	mpf_clear (term);
	mpf_clear (oterm);
	mpf_clear (acc);

	mpf_clear (one);
	mpf_clear (maxterm);
}

/* Implement Simon Plouffe odd-zeta sums */
static void fp_zeta_odd_helper (mpf_t zeta, unsigned int n, 
					 char *sdiv, char * spi, char * sminus, char * splus,  
					 unsigned int prec)
{
	mpf_t pi, pip, piterm, spos, sneg, spos_term, sneg_term, tmp;
	mpf_init (pi);
	mpf_init (pip);
	mpf_init (piterm);
	mpf_init (spos);
	mpf_init (sneg);
	mpf_init (spos_term);
	mpf_init (sneg_term);
	mpf_init (tmp);

	mpf_t div, c_pi, c_plus, c_minus;
	mpf_init (div);
	mpf_init (c_pi);
	mpf_init (c_plus);
	mpf_init (c_minus);
	
	mpf_set_str (div, sdiv, 10);
	mpf_set_str (c_pi, spi, 10);
	mpf_set_str (c_plus, splus, 10);
	mpf_set_str (c_minus, sminus, 10);
	
	fp_ess (spos, sneg, n, prec);
	mpf_mul (spos_term, spos, c_plus);
	mpf_mul (sneg_term, sneg, c_minus);
			  
	fp_pi (pi);
	mpf_pow_ui (pip, pi, n);
	mpf_mul (piterm, pip, c_pi);

	mpf_set (tmp, piterm);
	mpf_sub (zeta, tmp, spos_term);
	mpf_sub (tmp, zeta, sneg_term);
	mpf_div (zeta, tmp, div);
	
	mpf_clear (pi);
	mpf_clear (pip);
	mpf_clear (piterm);
	mpf_clear (spos);
	mpf_clear (sneg);
	mpf_clear (spos_term);
	mpf_clear (sneg_term);
	mpf_clear (tmp);
	
	mpf_clear (div);
	mpf_clear (c_pi);
	mpf_clear (c_plus);
	mpf_clear (c_minus);
}

int fp_zeta_odd_plouffe (mpf_t zeta, unsigned int n, unsigned int prec)
{
	int have_val = 1;
	switch (n)
	{
#ifdef PLOUFFE_ORIGINALS
		case 11: 
			fp_zeta_odd_helper (zeta, 11, "425675250", "1453", "851350500", "0", prec); 
			break;
		case 13: 
			fp_zeta_odd_helper (zeta, 13, "257432175", "89", "514926720", "62370", prec); 
			break;
		case 15: 
			fp_zeta_odd_helper (zeta, 15, "390769879500", "13687", "781539759000", "0", prec); 
			break;
		case 17: 
			fp_zeta_odd_helper (zeta, 17, "1904417007743250", "6758333", "3808863131673600", "29116187100", prec); 
			break;
		case 19: 
			fp_zeta_odd_helper (zeta, 19, "21438612514068750", "7708537", "42877225028137500", "0", prec); 
			break;
		case 21: 
			fp_zeta_odd_helper (zeta, 21, "1881063815762259253125", "68529640373", "3762129424572110592000", "1793047592085750", prec); 
			break;
#endif


		case 3:
			fp_zeta_odd_helper (zeta, 3, "180", "7", "360", "0", prec);
			break;

		case 5:
			fp_zeta_odd_helper (zeta, 5, "1470", "5", "3024", "84", prec);
			break;

		case 7:
			fp_zeta_odd_helper (zeta, 7, "56700", "19", "113400", "0", prec);
			break;

		case 9:
			fp_zeta_odd_helper (zeta, 9, "18523890", "625", "37122624", "74844", prec);
			break;

		case 11:
			fp_zeta_odd_helper (zeta, 11, "425675250", "1453", "851350500", "0", prec);
			break;

		case 13:
			fp_zeta_odd_helper (zeta, 13, "257432175", "89", "514926720", "62370", prec);
			break;

		case 15:
			fp_zeta_odd_helper (zeta, 15, "390769879500", "13687", "781539759000", "0", prec);
			break;

		case 17:
			fp_zeta_odd_helper (zeta, 17, "1904417007743250", "6758333", "3808863131673600", "29116187100", prec);
			break;

		case 19:
			fp_zeta_odd_helper (zeta, 19, "21438612514068750", "7708537", "42877225028137500", "0", prec);
			break;

		case 21:
			fp_zeta_odd_helper (zeta, 21, "1881063815762259253125", "68529640373", "3762129424572110592000", "1793047592085750", prec);
			break;

		case 23:
			fp_zeta_odd_helper (zeta, 23, "1211517431782539131250", "4472029801", "2423034863565078262500", "0", prec);
			break;

		case 25:
			fp_zeta_odd_helper (zeta, 25, "6948173623016040171631875", "2598638688071", "13896347660226074115072000", "414193993771808250", prec);
			break;

		case 27:
			fp_zeta_odd_helper (zeta, 27, "3952575621190533915703125", "149780635937", "7905151242381067831406250", "0", prec);
			break;

		case 29:
			fp_zeta_odd_helper (zeta, 29, "42344185423359347502790906715625", "162580897794660958", "84688371004458264623668408320000", "157739569618086594888750", prec);
			break;

		case 31:
			fp_zeta_odd_helper (zeta, 31, "28870481903812321637757079687500", "11231299844779783", "57740963807624643275514159375000", "0", prec);
			break;

		case 33:
			fp_zeta_odd_helper (zeta, 33, "17162190941764356274316709924901406250", "676470671886391879633", "34324381887524626998988066443264000000", "3995914450354646593461187500", prec);
			break;

		case 35:
			fp_zeta_odd_helper (zeta, 35, "923465669416292826066116829424218750", "3688053840923281541", "1846931338832585652132233658848437500", "0", prec);
			break;

		case 37:
			fp_zeta_odd_helper (zeta, 37, "3480645953760541547425811579090394140625", "1408434329374922032349", "6961291907571733063103925476843520000000", "50649968252302318662731718750", prec);
			break;

		case 39:
			fp_zeta_odd_helper (zeta, 39, "64875239172012679286579449799415644531250", "2659842854283579394387", "129750478344025358573158899598831289062500", "0", prec);
			break;

		case 41:
			fp_zeta_odd_helper (zeta, 41, "12967172230363787667401358845511389649052092451171875", "53866969189211783266383835533253", "25934344460739368914866833166704469676132761600000000", "11793580064115475681690378028576697656250", prec);
			break;

		case 43:
			fp_zeta_odd_helper (zeta, 43, "2919353325120984561556431951248804296083984375", "1228751826452728351300837", "5838706650241969123112863902497608592167968750", "0", prec);
			break;

		case 45:
			fp_zeta_odd_helper (zeta, 45, "25187657828037231081234525437683315511849888323792716796875", "1074151540472820600753617135934307286", "50375315656075893914882414558229863917282823887257600000000", "1431752413363682863232893583047239672166406250", prec);
			break;

		case 47:
			fp_zeta_odd_helper (zeta, 47, "15630294667467231804893395882010267487892367285156250", "67537532722660373286810600661", "31260589334934463609786791764020534975784734570312500", "0", prec);
			break;

		case 49:
			fp_zeta_odd_helper (zeta, 49, "125950123387606530332240169464377820361081282732372746373046875", "55141284330294633162607354950945193883", "251900246775213508129219880178582249920021784875591270400000000", "447464739541249826609197859219410845777653906250", prec);
			break;

		case 51:
			fp_zeta_odd_helper (zeta, 51, "5669518082718741943709352250640941481892180509033203125", "251492292317888012003479295207", "11339036165437483887418704501281882963784361018066406250", "0", prec);
			break;

		case 53:
			fp_zeta_odd_helper (zeta, 53, "19874174510194707355877113706035736161286459308343574747132284912109375", "89323943498389182315845947384336698100190998", "39748349020389419124707421860319918570808355970919336909209600000000000", "4412953194448248446248235437354232187414945030175781250", prec);
			break;

		case 55:
			fp_zeta_odd_helper (zeta, 55, "55921013802510257943421954165936644900510909734004056396484375", "25465609788816025420512226447159951", "111842027605020515886843908331873289801021819468008112792968750", "0", prec);
			break;

		case 57:
			fp_zeta_odd_helper (zeta, 57, "22954447164806465666694724159311794529795895568196223401126410316992032470703125", "1059122358196688900203789673076338001301584601329842", "45908894329612931651946395910854437757981681764417748097577037463552000000000000", "318556947592230848698389890628025301295324216829567935058593750", prec);
			break;

		case 59:
			fp_zeta_odd_helper (zeta, 59, "57522413794274203484482809918955361109315852011326949561410858154296875", "268916007610453025823381928132011055435166", "115044827588548406968965619837910722218631704022653899122821716308593750", "0", prec);
			break;

		case 61:
			fp_zeta_odd_helper (zeta, 61, "2250150739271701988086825366594937743084692866382599667359017830172943115234375", "1065838541236193393315346133195915243192115099012", "4500301478543403978125345388231641365295026250831015394540658360320000000000000", "1951694655041765879125640518065816059822622699974113769531250", prec);
			break;

		case 63:
			fp_zeta_odd_helper (zeta, 63, "463064280646029534081391924895270496216317118822260609653905317504882812500", "22223954766213317384532039590736747648635617", "926128561292059068162783849790540992432634237644521219307810635009765625000", "0", prec);
			break;

		case 65:
			fp_zeta_odd_helper (zeta, 65, "26404424191089874513907356219303296584761107108581828553012585401810261376879659509277343750", "128397633128226123041771286885908496599459788259186164309559", "52808848382179749029246099142333651014202872070357803507247465934466726253035520000000000000", "1431386703727057844680657853194146401222295130846203499276200981445312500", prec);
			break;

		case 67:
			fp_zeta_odd_helper (zeta, 67, "9549406932246083469671716326468923081604465571994768757702470085151062011718750", "4704971228496213648399974101471098623989701629", "19098813864492166939343432652937846163208931143989537515404940170302124023437500", "0", prec);
			break;

		case 69:
			fp_zeta_odd_helper (zeta, 69, "353445647207320312589501556954581303661421782728654161427946569538661284243177326949461956787109375", "17644260606276991066616771034325190251496670227214096141913859099", "706891294414640625180200634342120635660287073156637984171165512058662587799658225991680000000000000", "1197520432958028337443507699329661315272372981340019313303572092756086425781250", prec);
			break;

		case 71:
			fp_zeta_odd_helper (zeta, 71, "239800981547812029236551440284118221439891285487445040953972914548865012368257872009277343750", "1212919664600259084164537721476067892498197548805305939573", "479601963095624058473102880568236442879782570974890081907945829097730024736515744018554687500", "0", prec);
			break;

		case 73:
			fp_zeta_odd_helper (zeta, 73, "1291329911419567870960140556214255029201526083431909096561562593514463691075674732164924032352752685546875", "661787479183328801575691663514487261630136303487249188772092842108801", "2582659822839135741920554562173698885802190035618893569495071298019577925088504114397052928000000000000000", "273449745188827399137868755075376371946110990650542937154650067204863294494628906250", prec);
			break;

		case 75:
			fp_zeta_odd_helper (zeta, 75, "9223114674915847278328901549389162363072741749517116959768189021110192783394533538818359375", "478915836659382129612763840358992819232085424547810249", "18446229349831694556657803098778324726145483499034233919536378042220385566789067077636718750", "0", prec);
			break;

		case 77:
			fp_zeta_odd_helper (zeta, 77, "256915999653722002385420620581925893062631689304446702874814171671741700710513130338502685214062203340984344482421875", "1351677223440059667534026579871660256551221213825244302327592126797148471501838", "513831999307444004770844641418795271767626655863931280554152862827608390064812095216024486385863884800000000000000000", "3400254943485642363277255037874804524519484124988643785834539019115957739478118031311035156250", prec);
			break;

		case 79:
			fp_zeta_odd_helper (zeta, 79, "49730890667110062063780417096524621345713589999929696190237430116881950239082253054715538024902343750", "26509915083092912315730293342898546025433245995877421301815689", "99461781334220124127560834193049242691427179999859392380474860233763900478164506109431076049804687500", "0", prec);
			break;

		case 81:
			fp_zeta_odd_helper (zeta, 81, "10461016206657763952900402057947838221607973652032616773991058892666763307834520803510319865055904668370463053417205810546875", "565010122751068805311901539647748466645918598089367030993202855801018781817131385529", "20922032413315527905800812769045470197470074921363025449977418729263206306364731482648491918301061334433792000000000000000000", "8653149793754254127617297791901995300943929679690695689875627852188189251997692865893165588378906250", prec);
			break;

		case 83:
			fp_zeta_odd_helper (zeta, 83, "26003329595921675245471914390865080409321974652689816456520197013016076321011121153497045964145660400390625", "142302185198752633951003259184526109623109830283133673011603187937", "52006659191843350490943828781730160818643949305379632913040394026032152642022242306994091928291320800781250", "0", prec);
			break;

		case 85:
			fp_zeta_odd_helper (zeta, 85, "38887084245571232615762239155351229364466156226547279826038447112434066271185219471082222529642634439158386222743930816650390625", "21561960394860281920216321561911380255606766990072138345979703279525004543524413853342", "77774168491142465231524480321117594144504021436180079970640469144174461457175404223268649384341238425780748288000000000000000000", "2010415135415571708983085520318563574919306328914804965281104204325055969547463975842512138366699218750", prec);
			break;

		case 87:
			fp_zeta_odd_helper (zeta, 87, "5610873017648410876525990935061819676319817885052156089528252688404559699557514116011091518472362041473388671875", "315219849778284027953565106963657850893545087258111919846837707734741", "11221746035296821753051981870123639352639635770104312179056505376809119399115028232022183036944724082946777343750", "0", prec);
			break;

		case 89:
			fp_zeta_odd_helper (zeta, 89, "11983306424292523563272487428740775172961865538497429764189823547419372730118055752847146537708032471764200256195786163951258716583251953125", "68211882774187975535408148026046986746381545623428772099736410778808772736067950451251376306862", "23966612848585047126544974896201701707084680777660422979512776209763923457223152291551015253332510812672596352734421057536000000000000000000", "38720151361160949700665563451133129114925177996987040785856722177916445869144195840342848729633482566833496093750", prec);
			break;

		case 91:
			fp_zeta_odd_helper (zeta, 91, "98850114315536531663646110435038604330731448188564979377792850094858885927168976360594071271568450890002186298370361328125", "57011281443492086454066209812376204392756637590670916128169677803388899195966", "197700228631073063327292220870077208661462896377129958755585700189717771854337952721188142543136901780004372596740722656250", "0", prec);
			break;

		case 93:
			fp_zeta_odd_helper (zeta, 93, "204360232798243459477834020425322036897519411259299200372857852083551660110355253791974906854958627963075658651905457497836688137936115264892578125", "11942088084965907072485663838724465493618501820708239537797531067414422093362602012354356786198904884", "408720465596486918955668040891914294143556683991979354118834993008516690801243457625260357306325622591850127185391739613269196800000000000000000000", "41270220348517861473380953373119288841413370580532950041310543596408366665698809881580824617595820524127769470214843750", prec);
			break;

		case 95:
			fp_zeta_odd_helper (zeta, 95, "10413011211460717676250380323411563464940191729376966312797456539253869553468224695433458790140083700816192327770504951477050781250", "61653888415385795845515847826128391067032235441708936404643745680146263872104540729", "20826022422921435352500760646823126929880383458753932625594913078507739106936449390866917580280167401632384655541009902954101562500", "0", prec);
			break;

		case 97:
			fp_zeta_odd_helper (zeta, 97, "1026678332572275501552243051603595456179421058408388098350756910116148274900440004795140240146693296274946945963224938704627728474235150814056396484375", "615912172499623151056497876552021098489602355300407939720593786636702371108898393310875263926765893647", "2053356665144551003104486103220149414739702368026464258237535687690778991748288527895268450557922715410934530844309268265024094208000000000000000000000", "12958502380860251209688061536021867458482441947408518304987970264536122861040638917859390855768637259529698371887207031250", prec);
			break;

		case 99:
			fp_zeta_odd_helper (zeta, 99, "138242329869824150846594057174219062793343552968315549999019830815965948719447425713305868601176227514261236650225818157196044921875", "8402832178640067318790360222809790035833271338875779980075463442970968353567431741", "276484659739648301693188114348438125586687105936631099998039661631931897438894851426611737202352455028522473300451636314392089843750", "0", prec);
			break;

		case 101:
			fp_zeta_odd_helper (zeta, 101, "2419351059238966794264065234862224995465261436040175798246839347584984396457333341690077573976455955796273704240226369594280227140212861819492280483245849609375", "14899915005118230329522234211978534030191269202300034004489347623029060683618344857609303598783891763287037838", "4838702118477933588528130469726358522397161536334908646728357564776512561849762720980319611363558624681176873630030125885314483551941427200000000000000000000000", "1908531466638664254557050234678869606543768935096037600164463410646713088629465149577386696754029271515703561015439033508300781250", prec);
			break;

		case 103:
			fp_zeta_odd_helper (zeta, 103, "8197924577963036738779523235993054026338412855769777723411224872537802193027951967573559770704977805441824763159729318960607051849365234375", "5115511590011241057447329371766791619756062974358379446868597834892199267027172693354161", "16395849155926073477559046471986108052676825711539555446822449745075604386055903935147119541409955610883649526319458637921214103698730468750", "0", prec);
			break;

		case 105:
			fp_zeta_odd_helper (zeta, 105, "33476488358106009825763023779631159333930662750380792623988169928684683541846583233039949042510088277372839071980087962351119169683052757572819981952352344989776611328125", "2116534265787916485294070652838082080793486125598617261243576592838819586359680190232924753656090810968893016676643182", "66952976716212019651526047559263969186168187761577827092747575098869588476838198884290140630310099544251049313116201021657262466921372241384334950400000000000000000000000", "1650518306862260816241844771235241500221393145032418210242545289922989505371169156025096955024127555266726238694986495295310020446777343750", prec);
			break;

		case 107:
			fp_zeta_odd_helper (zeta, 107, "2783058325288414036849712398236214583698258217456119644924215406100220131744523364320583550775783967400316521027570187455100892106179893016815185546875", "17828219473673421907944196084848100601631041343219326565113363472678592211417261132053010397435946", "5566116650576828073699424796472429167396516434912239289848430812200440263489046728641167101551567934800633042055140374910201784212359786033630371093750", "0", prec);
			break;

		case 109:
			fp_zeta_odd_helper (zeta, 109, "2770331665932256818162743030592561365234229703112462320967361126017613656158622946584043403933481618972813773164102032875721296469999862273549059875039696490764617919921875", "1798115913378776908280510846581495927392919500045134728012557199809377541708858010727893155140611067717369930678676476", "5540663331864513636325486061185131267211997491200411704593289988764979140338855715882371890888045154659787061831237323833554579690290949200850095964160000000000000000000000", "8536743538084975487062658567736729751828021609822714285083021081916714159515503033258082111986750291224653751976214080607018470764160156250", prec);
			break;

		case 111:
			fp_zeta_odd_helper (zeta, 111, "14154722304207990055970853205231108472106529924446704914669336044188681858792198003618005319551095388809433241618523006647676861851113894276320934295654296875", "930866768598038555037745721541506786862470847228989799903334802511940872169956683902397476952340874379", "28309444608415980111941706410462216944213059848893409829338672088377363717584396007236010639102190777618866483237046013295353723702227788552641868591308593750", "0", prec);
			break;

		case 113:
			fp_zeta_odd_helper (zeta, 113, "535988016715013803215885177694905048863337478198426795760058206015602811630153582231911404971143236251176709142305551185761870020006558017553544133682564617258139066113531589508056640625", "3571424780838804649141752172968771217680454853370454729413198784149643035808789811076069610075146232317965831107065035516936474642", "1071976033430027606431770355389810200954212078813991794526733927062419132270996124333976445113230033066717508493692197433458049629408308677073541512606489446973440000000000000000000000000", "103227537122417138203006617515031213509010688959870153635170943560564364090209081095061934309589395192641966453245241360212457161867772936820983886718750", prec);
			break;

		case 115:
			fp_zeta_odd_helper (zeta, 115, "3071172780353379457032523843415732210823162649955964877499867887705336171175304348586221450337570021206447556909145132762843916595594523182667791843414306640625", "2073437420382515647798623037303990230623011183648874780288500296064825521017059758527851353888748509478", "6142345560706758914065047686831464421646325299911929754999735775410672342350608697172442900675140042412895113818290265525687833191189046365335583686828613281250", "0", prec);
			break;

		case 117:
			fp_zeta_odd_helper (zeta, 117, "21498736624687227040086513735664515356621082823324464377234984302189136308352789944695014299533816480962961428270060655639882451548138069036834366995973401718165905131535342909395694732666015625", "1470617995945929340823900939912466478922162017902177245157936266684462755304379703178802123827766074689212758305303974887068597968722324", "42997473249374454080173027471329030972023794596522064601278273085897596918512021096333316097845847228622530016105717753896003333477478893483777645464708595559741731635200000000000000000000000000", "258781628949873135846808304481519324301806441206943287498778214266696607159565596442616238430381202755410108911472761792123409921372129314181208610534667968750", prec);
			break;

		case 119:
			fp_zeta_odd_helper (zeta, 119, "181961744064274353998732174876458083192401580310810970658317353559398594886909058323836349301639835028493133397370199672540329374986916978202036009397671930491924285888671875", "1261151562313641703816834464893880014977525934624571730625687371825311625760783202572888567761894854518812674014678", "363923488128548707997464349752916166384803160621621941316634707118797189773818116647672698603279670056986266794740399345080658749973833956404072018795343860983848571777343750", "0", prec);
			break;
			
		default:
			have_val = 0;
	}
	return have_val;
}

/* ======================================================================= */
/* Brute force summation of zeta values */

void fp_zeta_brute (mpf_t zeta, unsigned int s, int prec)
{
	unsigned long int us = s;

	/* Set up cache of what we've computed so far */
	static int cache_size = 0;
	static mpf_t *zeta_cache = NULL;
	static int *zprec= NULL;
	static unsigned int *last_term= NULL;

	if (s >= cache_size)
	{
		int newsize = (3*s)/2+20;
		zeta_cache = (mpf_t *) realloc (zeta_cache, newsize*sizeof (mpf_t));
		zprec = (int *) realloc (zprec, newsize*sizeof (int));
		last_term = (unsigned int *) realloc (last_term, newsize*sizeof (unsigned int));
		
		int i;
		for (i=cache_size; i<newsize; i++)
		{
			zprec[i] = 0;
		}
		cache_size = newsize;
	}

	if (s<2) return;
	
	/* Lets see if we can get lucky with the cache. */
	if (prec < zprec[s])
	{
		mpf_set (zeta, zeta_cache[s]);
		return;
	}

	/* Initialize the cache line, if needed */
	if (0 == zprec[s])
	{
		mpf_init (zeta_cache[s]);
		mpf_set_ui (zeta_cache[s], 1);
		last_term[s] = 2;
		zprec[s] = prec;
	}

	/* If we are here, well have to compute values using brute force */
	mpf_t acc;
	mpf_t term;
	mpf_t en;
	mpf_t inv;
	
	mpf_init (acc);
	mpf_init (term);
	mpf_init (en);
	mpf_init (inv);
	
	mpf_set_ui (zeta, 1);

	/* Compute number of terms to be carried out.
	 * If we want error t be less than epsilon, 
	 * then must sum to epsilon=N^{1-s} or
	 * N = exp (ln(espilon) / (1-s))
	 * But epsilon = 10^{-prec} so
	 * Nmax = 10^{prec/(s-1)}
	 *
	 * Note that this provides a precise upper bound
	 * on the error term, for an s>7.
	 */
	double fprec = prec;
	fprec /= (double) (s-1);
	double fnmax = pow (10.0, fprec);
	if (1.0e9 < fnmax)
	{
		fprintf (stderr, "Sorry bucko, can't do it, you asked for zeta(%d) in %g digits\n", s, fnmax);
		return;
	}
	int nmax = (int) (fnmax+1.0);
	printf ("zeta(%d) to precision %d will require %d terms\n", s, prec, nmax);
	
	/* Start computations where we last left off. */
	mpf_set (zeta, zeta_cache[s]);
	int nstart = last_term[s];
	
	int n;
	for (n=nstart; n< nmax; n++)
	{
		mpf_set_ui (en, n);
		mpf_ui_div (inv, 1, en);  /* inv = 1/n */
		mpf_pow_ui (term, inv, us); /* term = 1/n^s */
		mpf_add (acc, zeta, term);
		mpf_set (zeta, acc);
	}

	/* cache the results */
	mpf_set (zeta_cache[s], zeta);
	last_term[s] = nmax;
	
	mpf_clear (acc);
	mpf_clear (term);
	mpf_clear (en);
	mpf_clear (inv);
}

/* ======================================================================= */
/* the d_k from the borwein 1995 paper */

static void fp_borwein_tchebysheff (mpf_t d_k, int n, int k)
{
	DECLARE_FP_CACHE (cache);

	/* cache the likeliest case: same value of n every time. */
	static int ncache = 0;
	if (n != ncache)
	{
		fp_one_d_cache_clear(&cache);
		ncache = n;
	}
	if ((0 == k) || (0 == n))
	{
		mpf_set_ui (d_k, 1); 
		return;
	}
	int hit = fp_one_d_cache_check (&cache, k);
	if (hit)
	{
		fp_one_d_cache_fetch (&cache, d_k, k);
		return;
	}

	mpz_t ifact;
	mpz_init (ifact);

	mpf_t term, fact, four;
	mpf_init (term);
	mpf_init (fact);
	mpf_init (four);

	mpf_set_ui (d_k, 0);
	mpf_set_ui (four, 1);
	int i;
	for (i=0; i<=n; i++)
	{
		i_factorial (ifact, n+i-1);
		mpf_set_z (term, ifact);
		i_factorial (ifact, n-i);
		mpf_set_z (fact, ifact);
		mpf_div (term, term, fact);
		i_factorial (ifact, 2*i);
		mpf_set_z (fact, ifact);
		mpf_div (term, term, fact);
		mpf_mul (term, term, four);
		mpf_mul_ui(term, term, n);

		mpf_add (d_k, d_k, term);

		fp_one_d_cache_store (&cache, d_k, i, 1);

		mpf_mul_ui (four, four, 4);
	}

	mpf_clear (fact);
	mpf_clear (term);
	mpf_clear (four);
	mpz_clear (ifact);

	fp_one_d_cache_fetch (&cache, d_k, k);
}

void fp_borwein_zeta (mpf_t zeta, unsigned int s, int prec)
{
	double nterms = 0.69 + 2.302585093 * prec;
	// Huh? whazzup with the gamma ??
	// nterms -=  s * log(s) -s;
	nterms *= 0.567296329;
	int n = nterms+1.0;

	mpz_t ip;
	mpz_init (ip);

	mpf_t d_n, po, term, twon;
	mpf_init (d_n);
	mpf_init (po);
	mpf_init (term);
	mpf_init (twon);

	fp_borwein_tchebysheff (d_n, n, n);

	mpf_set_ui (zeta, 0);
	int k;
	for (k=0; k<n; k++)
	{
		fp_borwein_tchebysheff (term, n, k);
		mpf_sub (term, term, d_n); 

		// i_pow (ip, k+1, s);
		mpz_ui_pow_ui (ip, k+1, s);
		mpf_set_z (po, ip);
		mpf_div (term, term, po);

		if (k%2)
		{
			mpf_sub(zeta, zeta, term);
		}
		else
		{
			mpf_add(zeta, zeta, term);
		}
	}
	mpf_div (zeta, zeta, d_n);
	mpf_neg (zeta, zeta);

	mpf_set_ui (twon, 1);
	mpf_div_2exp (twon, twon, s-1);
	
	mpf_set_ui (term, 1);
	mpf_sub (term, term, twon);
	
	mpf_div (zeta, zeta, term);

	mpz_clear (ip);
	mpf_clear (twon);
	mpf_clear (term);
	mpf_clear (po);
	mpf_clear (d_n);
}

/* ======================================================================= */
/* fp_zeta
 * Floating-point-valued Riemann zeta for positive integer arguments 
 * return value placed in the arg "zeta".
 *
 * Carries out the math to "prec" decimal digits. Uses a combined
 * algorithm: 
 * 
 * For even "n", computes an "exact" result be using recursion
 * to get the Bernoulli numbers, and then working off of those.
 *
 * For odd "n", this uses a fast convergent sum based on consts 
 * from Simon Plouffe for low odd values of "n" (n less than 120). 
 *
 * For large odd "n", performs the brute-force sum. This works,
 * and works quite well, since the sum converges decently when "n"
 * is large.
 *
 */
void fp_zeta (mpf_t zeta, unsigned int s, int prec)
{
	DECLARE_FP_CACHE (cache);
	if (2>s)
	{
		fprintf (stderr, "Domain error, asked for zeta(%d)\n", s);
		mpf_set_ui (zeta, 0);
		return;
	}

	int have_prec = fp_one_d_cache_check (&cache, s);
	if (have_prec >= prec)
	{
		fp_one_d_cache_fetch (&cache, zeta, s);
		return;
	}
	
	/* We've got exact results for even numbers */
	if (0 == s%2)
	{
		fp_zeta_even (zeta, s);
		fp_one_d_cache_store (&cache, zeta, s, 10111222);
		return;
	}
	
	/* Bump precision so as to increase cache hits on next go-around. */
	// prec += 100;

	int plo = 0;
	// plo = fp_zeta_odd_plouffe (zeta, s, prec);

	if (0 == plo)
	{
		// fp_zeta_brute (zeta, s, prec);
		fp_borwein_zeta (zeta, s, prec);
	}

	/* Save computed value to the cache. */
	fp_one_d_cache_store (&cache, zeta, s, prec);
}

/* ======================================================================= */
/* rough count of number of digits in a number */

static inline unsigned int num_digits (mpz_t num, mpz_t tmpa, mpz_t tmpb)
{
	unsigned int n=0;
	
	mpz_set (tmpb, num);
	while (1)
	{
		mpz_fdiv_q_ui (tmpa, tmpb, 100);
		mpz_set (tmpb, tmpa);
		if (0 == mpz_sgn  (tmpa)) break;
		n += 2;
	}
	return n;
}

/* ======================================================================= */
/* 
 * Compute a_sub_n
 * the w argument is for the power bit -- 
 */
void a_sub_n (mpf_t a_n, mpf_t w, unsigned int n, unsigned int prec)
{
	int k;
	mpf_t fbin, term, zt, ok, one, acc, zeta;
	mpf_t gam, wneg, wn;

	mpf_init (term);
	mpf_init (acc);
	mpf_init (zeta);
	mpf_init (zt);
	mpf_init (ok);
	mpf_init (one);
	mpf_init (fbin);
	mpf_init (gam);
	mpf_init (wneg);
	mpf_init (wn);
	
	mpz_t tmpa, tmpb;
	mpz_init (tmpa);
	mpz_init (tmpb);
	
	mpf_set_ui (one, 1);

	mpz_t ibin;
	mpz_init (ibin);
	mpf_set_ui (a_n, 0);

	mpf_neg (wneg, w);
	mpf_set (wn, wneg);

	int maxbump = 0;
	for (k=1; k<= n; k++)
	{
		/* Commpute the binomial */
		i_binomial (ibin, n, k);
		mpf_set_z (fbin, ibin);

		/* The terms will have alternating signs, and
		 * will mostly cancel one-another. Thus, we need 
		 * to increase precision for those terms with the 
		 * largest binomial coefficients. This is will
		 * increase precision for the killer terms, 
		 * while keeping the others in bearable range,
		 * in terms to cpu time consumed.
		 */
		int ndigits = num_digits (ibin, tmpa,tmpb);
		if (maxbump < ndigits) maxbump = ndigits;

		/* compute 1/k - zeta (k+1)/(k+1) */
		fp_zeta (zeta, k+1, prec+ndigits);
		// fp_hasse_zeta (zeta, k+1, prec+ndigits);

		mpf_div_ui (zt, zeta, k+1);
		mpf_div_ui (ok, one, k);
		mpf_sub (term, ok, zt);

		mpf_mul (zeta, term, fbin);

#define W_IS_EQUAL_TO_ONE 1
#if W_IS_EQUAL_TO_ONE
		if (k%2) mpf_neg (term, zeta);
		else mpf_set (term, zeta);
#else 
		mpf_mul (term, wn, zeta);
		mpf_mul (zt, wn, wneg);
		mpf_set (wn, zt);
#endif
		
		mpf_add (acc, a_n, term);
		mpf_set (a_n, acc);
	}

	/* add const terms */
	mpf_add_ui (term, a_n, 1);
	fp_euler_mascheroni (gam);
	mpf_sub (a_n, term, gam);

	/* subtract 1/2(n+1) */
	mpf_div_ui (ok, one, 2*(n+1));
#if W_IS_EQUAL_TO_ONE
#else
	mpf_div (zt, ok, w);
	mpf_set (ok, zt);
#endif
	mpf_sub (term, a_n, ok);
	mpf_set (a_n, term);
	
	mpf_clear (term);
	mpf_clear (acc);
	mpf_clear (zeta);
	mpf_clear (zt);
	mpf_clear (ok);
	mpf_clear (one);
	mpf_clear (fbin);
	mpf_clear (gam);
	mpf_clear (wneg);
	mpf_clear (wn);

	mpz_clear (ibin);
	mpz_clear (tmpa);
	mpz_clear (tmpb);

	// printf ("# max precision bump=%d\n", maxbump);
}

void a_bound_n (mpf_t b_n, unsigned int n)
{
	mpf_t en, sq_en;
	mpf_init (en);
	mpf_init (sq_en);

	mpf_set_ui (en, n+1);
	mpf_sqrt (sq_en, en);
	mpf_neg (en, sq_en);

	mpf_clear (en);
	mpf_clear (sq_en);
}

/* ======================================================================= */
/* 
 * Compute b_sub_n
 */
void b_sub_n (mpf_t b_n, unsigned int n, unsigned int prec)
{
	DECLARE_FP_CACHE (cache);
	if (0 == n)
	{
		mpf_set_ui (b_n, 1);
		mpf_div_ui (b_n, b_n, 2);
		return;
	}

	int have_prec = fp_one_d_cache_check (&cache, n);
	if (have_prec >= prec)
	{
		fp_one_d_cache_fetch (&cache, b_n, n);
		return;
	}

	if (1 == n)
	{
		mpf_set_ui (b_n, 1);
		mpf_div_ui (b_n, b_n, 2);
		mpf_t gam;
		mpf_init (gam);
		fp_euler_mascheroni (gam);
		mpf_sub(b_n, b_n, gam);
		mpf_clear (gam);
		fp_one_d_cache_store (&cache, b_n, n, 10111222);
		return;
	}
	
	mpz_t ibin;
	mpz_init (ibin);
	
	mpf_t bin, zeta;
	mpf_init (bin);
	mpf_init (zeta);

	mpf_set_si (b_n, -1);
	mpf_div_ui (b_n, b_n, 2);

	int k;
	for (k=2; k<=n; k++)
	{
		i_binomial (ibin, n, k);
		mpf_set_z (bin, ibin);
		fp_zeta (zeta, k, prec);
		mpf_mul (zeta, zeta, bin);
		if (k%2)
		{
			mpf_sub (b_n, b_n, zeta);
		}
		else
		{
			mpf_add (b_n, b_n, zeta);
		}
	}

	/* now for the oddball terms */
	fp_euler_mascheroni (bin);
	mpf_set_ui (zeta, 1);
	mpf_sub (zeta, zeta, bin);
	fp_harmonic (bin, n-1);
	mpf_sub (zeta, zeta, bin);
	mpf_mul_ui (zeta, zeta, n);
	mpf_add (b_n, b_n, zeta);

	mpf_clear (bin);
	mpf_clear (zeta);
	mpz_clear (ibin);

	fp_one_d_cache_store (&cache, b_n, n, prec);
}

/* ======================================================================= */
/* compute a_sub_s for complex-valued s
 */
void a_sub_s (mpf_t re_a, mpf_t im_a, double re_s, double im_s, unsigned int prec)
{
	int k;
	mpf_t rebin, imbin, term, zt, ok, one, racc, iacc, rzeta, izeta;
	mpf_t gam;

	mpf_init (term);
	mpf_init (racc);
	mpf_init (iacc);
	mpf_init (rzeta);
	mpf_init (izeta);
	mpf_init (ok);
	mpf_init (one);
	mpf_init (zt);
	mpf_init (rebin);
	mpf_init (imbin);
	mpf_init (gam);
	
	mpz_t tmpa, tmpb;
	mpz_init (tmpa);
	mpz_init (tmpb);
	
	mpf_set_ui (one, 1);
	mpf_set_ui (re_a, 0);
	mpf_set_ui (im_a, 0);

	int n = 1500;  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	for (k=1; k<= n; k++)
	{
		/* Commpute the binomial */
		c_binomial (rebin, imbin, re_s, im_s, k);

		/* compute 1/k - zeta (k+1)/(k+1) */
		int ndigits = 0;
		fp_zeta (rzeta, k+1, prec+ndigits);
		mpf_div_ui (zt, rzeta, k+1);
		mpf_div_ui (ok, one, k);
		mpf_sub (term, ok, zt);

		mpf_mul (rzeta, term, rebin);
		mpf_mul (izeta, term, imbin);

		if (k%2)
		{ 
			mpf_sub (racc, re_a, rzeta);
			mpf_sub (iacc, im_a, izeta);
		}
		else 
		{
			mpf_add (racc, re_a, rzeta);
			mpf_add (iacc, im_a, izeta);
		}
		
		mpf_set (re_a, racc);
		mpf_set (im_a, iacc);
	}

	/* add const terms */
	mpf_add_ui (term, re_a, 1);
	fp_euler_mascheroni (gam);
	mpf_sub (re_a, term, gam);

	/* subtract 1/2(s+1) */
	double rex = 2.0*(re_s +1.0);
	double imx = 2.0*im_s;
	double den = rex*rex + imx*imx;
	rex = rex / den;
	imx = -imx / den;

	mpf_set_d (ok, rex);
	mpf_sub (term, re_a, ok);
	mpf_set (re_a, term);
	
	mpf_set_d (ok, imx);
	mpf_sub (term, im_a, ok);
	mpf_set (im_a, term);
	
	mpf_clear (term);
	mpf_clear (racc);
	mpf_clear (iacc);
	mpf_clear (rzeta);
	mpf_clear (izeta);
	mpf_clear (ok);
	mpf_clear (one);
	mpf_clear (zt);
	mpf_clear (rebin);
	mpf_clear (imbin);
	mpf_clear (gam);

	mpz_clear (tmpa);
	mpz_clear (tmpb);

}

/* ======================================================================= */
/* compute b_sub_s for complex-valued s
 */
void b_sub_s (mpf_t re_b, mpf_t im_b, double re_s, double im_s, unsigned int prec, int nterms)
{
	int k;
	mpf_t rebin, imbin, term, ok, one, racc, iacc, rzeta, izeta;
	mpf_t gam;

	mpf_init (term);
	mpf_init (racc);
	mpf_init (iacc);
	mpf_init (rzeta);
	mpf_init (izeta);
	mpf_init (ok);
	mpf_init (one);
	mpf_init (rebin);
	mpf_init (imbin);
	mpf_init (gam);
	
	mpf_set_ui (one, 1);
	mpf_set_ui (re_b, 0);
	mpf_set_ui (im_b, 0);
	fp_euler_mascheroni (gam);

	int n = 650;  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	n = nterms;
	for (k=2; k<= n; k++)
	{
		/* Commpute the binomial */
		c_binomial (rebin, imbin, re_s, im_s, k);

// printf ("duude s= (%g %g) k=%d bin=(%g %g)\n", re_s, im_s, k, mpf_get_d(rebin), mpf_get_d(imbin));

		/* compute zeta (k)- 1/(k-1) - gamma */
		int ndigits = 0;
		fp_zeta (rzeta, k, prec+ndigits);
		mpf_div_ui (ok, one, k-1);
		mpf_sub (term, rzeta, ok);
		mpf_sub (term, term, gam);

		mpf_mul (rzeta, term, rebin);
		mpf_mul (izeta, term, imbin);

		if (k%2)
		{ 
			mpf_sub (racc, re_b, rzeta);
			mpf_sub (iacc, im_b, izeta);
		}
		else 
		{
			mpf_add (racc, re_b, rzeta);
			mpf_add (iacc, im_b, izeta);
		}
		
		mpf_set (re_b, racc);
		mpf_set (im_b, iacc);
	}

	/* add const terms */
	mpf_sub (re_b, re_b, gam);

	/* add 1/2 */
	mpf_div_ui (ok, one, 2);
	mpf_add (re_b, re_b, ok);
	
	mpf_clear (term);
	mpf_clear (racc);
	mpf_clear (iacc);
	mpf_clear (rzeta);
	mpf_clear (izeta);
	mpf_clear (ok);
	mpf_clear (one);
	mpf_clear (rebin);
	mpf_clear (imbin);
	mpf_clear (gam);
}

/* ==================================================================== */
/* Return the Steiltjes constants */

#if THIS_WORKS_BUT_HAS_PRECISION_PROBLEMS
void stieltjes_gamma (mpf_t gam, int n)
{
	int k;

	mpz_t isb;
	mpz_init (isb);

	mpf_t term, sb;
	mpf_init (term);
	mpf_init (sb);

	mpf_set_ui (gam, 0);
	// XXXX precision violation !!
	for (k=n; k<n+260; k++)
	{
		b_sub_n (term, k, 460);
		i_stirbin_sum (isb, k,n);
		mpf_set_z (sb, isb);
		mpf_mul (term, term, sb);

		i_factorial (isb, k);
		mpf_set_z (sb, isb);
		mpf_div (term, term, sb);
		mpf_add (gam, gam, term);
	}
	i_factorial (isb, n);
	mpf_set_z (sb, isb);
	mpf_mul (gam, gam, sb);
	if (n%2) mpf_neg (gam, gam);

	mpf_clear (term);
	mpf_clear (sb);
	mpz_clear (isb);
}
#endif

/* ==================================================================== */
#ifdef TEST

main (int argc, char * argv[])
{
	char str[4000];

#ifdef FACT_TEST
	mpz_t fact;
	mpz_init (fact);

	i_factorial (fact, 5);
	mpz_get_str (str, 10, fact);
	printf ("fact = %s\n", str);
#endif

#ifdef I_BINOMIAL_TEST
	int n, k;
	mpz_t bin;
	mpz_init (bin);

	for (n=1; n<7; n++)
	{
		for (k=0; k<=n; k++)
		{
			i_binomial (bin, n ,k);
			mpz_get_str (str, 10, bin);
			printf ("bin (%d %d) = %s\n", n, k, str);
		}
		printf ("---\n");
	}
#endif

// #define I_STIRLING_TEST
#ifdef I_STIRLING_TEST
	int n, k;
	mpz_t sitrly;
	mpz_init (sitrly);

	for (n=0; n<21; n++)
	{
		for (k=0; k<=n; k++)
		{
			i_stirling_first (sitrly, n ,k);
			mpz_get_str (str, 10, sitrly);
			printf ("sitrly (%d %d) = %s\n", n, k, str);
		}
		printf ("---\n");
	}
#endif

// #define F_BINOMIAL_TEST
#ifdef F_BINOMIAL_TEST
	int n, k;
	mpf_t bin;
	mpf_init (bin);

	for (n=1; n<7; n++)
	{
		for (k=0; k<=n; k++)
		{
			fp_binomial (bin, (double)n ,k);
			printf ("bin (%d %d) = ", n, k);
			mpf_out_str (stdout, 10, 60, bin);
			printf ("\n");
		}
		printf ("---\n");
	}
#endif
	
	if (argc < 3)
	{
		fprintf (stderr, "Usage: %s [ndigits] [nterms]\n", argv[0]); 
		exit (1);
	}

	/* the decimal precison (number of decimal places) */
	int prec = atoi (argv[1]);

	/* number of an's to compute */
	int nterms = atoi (argv[2]);

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
	/* XXX this should be stirling's approx for binomial */ 
	int bits = (int) (v + 300 + 3*nterms);
	
	/* set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
#ifdef ZETA_STUFF
	mpf_t zeta;
	mpf_init (zeta);
	
	printf ("           000000000011111111112222222222333333333344444444445555555555666666666677777777778\n");
	printf ("           012345678901234567890123456789012345678901234567890123456789012345678901234567890\n");
	fp_zeta (zeta, 2, 13);
	fp_prt ("13 digs= ", zeta);
	fp_zeta (zeta, 8, 30);
	fp_prt ("30 digs= ", zeta);
	fp_zeta (zeta, 8, 40);
	fp_prt ("40 digs= ", zeta);
	fp_zeta (zeta, 8, 50);
	fp_prt ("50 digs= ", zeta);
	fp_zeta (zeta, 8, 60);
	fp_prt ("60 digs= ", zeta);
	fp_zeta (zeta, 8, 70);
	fp_prt ("70 digs= ", zeta);
	fp_zeta (zeta, 8, 80);
	fp_prt ("0 digs= ", zeta);
#endif
	
// #define TEST_ZETA
#ifdef TEST_ZETA
	mpf_t zeta;
	mpf_init (zeta);
	// fp_zeta_odd (zeta, 3, 180, 7, 360, 0, 60); 
	// fp_prt ("duude zeta3= ", zeta);
	int i;
	int pr = prec;
	for (i=3; i<nterms; i++ ) {
		fp_zeta (zeta, i, pr);
		printf ("char * zeta_%d_%d = \"", i, pr);
		mpf_out_str (stdout, 10, pr, zeta);
		printf ("\";\n");
		fflush (stdout);
	}
#endif

#ifdef TEST_EXP
	mpf_t ex, one;
	mpf_init (ex);
	mpf_init (one);
	mpf_set_ui(one, 1);
	fp_exp (ex, one, 50);
	fp_prt ("e= ", ex);
	mpf_clear (ex);
	mpf_clear(one);
#endif
	
#ifdef TEST_DIGIT_COUNT
	mpz_t ival, tmpa, tmpb;
	mpz_init (ival);
	mpz_init (tmpa);
	mpz_init (tmpb);
	mpz_set_ui (ival, 3000000);
	int nd = num_digits (ival, tmpa, tmpb);
	printf ("found %d digits\n", nd);
#endif

// #define TEST_BERNOULLI
#ifdef TEST_BERNOULLI
	mpq_t bern;
	mpq_init (bern);
	int n = 4;
	for (n=8; n<30; n++) 
	{
		q_bernoulli (bern, n);
		printf ("bernoulli (%d)= ", n);
		mpq_out_str (stdout, 10, bern);
		printf ("\n");
	}
#endif /* TEST_BERNOULLI */

#define TEST_STIELTJES
#ifdef TEST_STIELTJES
	mpf_t stie;
	mpf_init (stie);
	int i;
	for (i=0; i<40; i++ ) {
		stieltjes_gamma (stie, i);
		printf ("gamma[%d] = ", i);
		mpf_out_str (stdout, 10, 60, stie);
		printf (";\n");
		fflush (stdout);
	}
#endif

}
#endif

/* =============================== END OF FILE =========================== */


/*
 * mp_zeta.c
 *
 * High-precison Riemann zeta function, using the 
 * Gnu Multiple-precision library.
 *
 * Actually, high-precision a_s on the complex plane
 * XXXX actually, this is not the original source, use the 
 * other file mp_zeta.c  for teh more curently maintained 
 * version of these routines.  these are cut-n-pasted from 
 * there,for convenience
 * 
 * Linas Vepstas July 2005
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "brat.h"

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
	static i_cache name = {0, NULL, NULL, 0}

/** i_one_d_cache_check() -- check if mpz_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int i_one_d_cache_check (i_cache *c, unsigned int n)
{
	if (c->disabled) return 0;
	if ((n > c->nmax) || 0==n )
	{
		unsigned int newsize = (int) (1.5*n+1.0);
		c->cache = (mpz_t *) realloc (c->cache, newsize * sizeof (mpz_t));
		c->ticky = (char *) realloc (c->ticky, newsize * sizeof (char));

		unsigned int en;
		unsigned int nstart = c->nmax+1;
		if (0 == c->nmax) nstart = 0;
		for (en=nstart; en <newsize; en++)
		{
			mpz_init (c->cache[en]);
			c->ticky[en] = 0;
		}
		c->nmax = newsize-1;
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
/* Cache management */
/* pure cut-n-paste of he integer variant */

typedef struct {
	unsigned int nmax;
	mpq_t *cache;
	char *ticky;
} q_cache;

#define DECLARE_Q_CACHE(name)         \
	static q_cache name = {0, NULL, NULL}

/** q_one_d_cache_check() -- check if mpq_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int q_one_d_cache_check (q_cache *c, unsigned int n)
{
	if ((n > c->nmax) || 0==n )
	{
		unsigned int newsize = (int) (1.5*n+1.0);
		c->cache = (mpq_t *) realloc (c->cache, newsize * sizeof (mpq_t));
		c->ticky = (char *) realloc (c->ticky, newsize * sizeof (char));

		unsigned int en;
		unsigned int nstart = c->nmax+1;
		if (0 == c->nmax) nstart = 0;
		for (en=nstart; en <newsize; en++)
		{
			mpq_init (c->cache[en]);
			c->ticky[en] = 0;
		}
		c->nmax = newsize-1;
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
	static fp_cache name = {0, NULL, NULL}

/** fp_one_d_cache_check() -- check if mpf_t value is in the cache
 *  If there is a cached value, this returns the precision of the 
 *  value in the cache; else it returns zero.
 *  This assumes a 1-dimensional cache layout (simple array)
 */
int fp_one_d_cache_check (fp_cache *c, unsigned int n)
{
	if ((n > c->nmax) || 0==n )
	{
		unsigned int newsize = (int) (1.5*n+1.0);
		c->cache = (mpf_t *) realloc (c->cache, newsize * sizeof (mpf_t));
		c->precision = (int *) realloc (c->precision, newsize * sizeof (int));

		unsigned int en;
		unsigned int nstart = c->nmax+1;
		if (0 == c->nmax) nstart = 0;
		for (en=nstart; en <newsize; en++)
		{
			mpf_init (c->cache[en]);
			c->precision[en] = 0;
		}
		c->nmax = newsize-1;
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
	unsigned int i;
	for (i=0; i<c->nmax; i++)
	{
		c->precision[i] = 0;
	}
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
#define USE_LOCAL_FACTORIAL
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
#define i_binomial mpz_bin_uiui

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
	
	unsigned int istart = n-1;
	hit = fp_one_d_cache_check (&cache, istart);
	while (0 == hit && istart>1)
	{
		istart--;
		hit = fp_one_d_cache_check (&cache, istart);
	}

	unsigned int i;
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

/* ======================================================================= */
/**
 * fp_zeta_even - return the zeta value for even "n".
 * @prec - decimal places of precision to work to.
 *
 * Uses a fast algorithm to compute the zeta function for any 
 * even value of n. This is obtained by recursievly computing
 * the Bernoulli numbers, and multiplying by an appropriate factor
 * of pi and factorial. 
 */
void fp_zeta_even (mpf_t zeta, unsigned int n, int prec)
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
	
	fp_pi (pi, prec);
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
	fp_one_d_cache_check (&cache, n);
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
	int n = (int) (nterms+1.0);

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
		fp_zeta_even (zeta, s, prec);
		fp_one_d_cache_store (&cache, zeta, s, prec);
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
	fp_euler_mascheroni (gam,prec);
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
	fp_euler_mascheroni (gam, prec);

	int n = 650;  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	n = nterms;
	int downer  = 0;
	for (k=2; k<= n; k++)
	{
		/* Commpute the binomial */
		c_binomial (rebin, imbin, re_s, im_s, k);

// printf ("duude s= (%g %g) k=%d bin=(%g %g)\n", re_s, im_s, k, mpf_get_d(rebin), mpf_get_d(imbin));

		/* compute zeta (k)- 1/(k-1) */
		// fp_zeta_minus_pole (rzeta, k, prec);
		fp_zeta (rzeta, k, prec);
		mpf_div_ui (ok, one, k-1);
		mpf_sub (term, rzeta, ok);

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
#if 1
		double rt = mpf_get_d (rzeta);
		double it = mpf_get_d (izeta);
		double ra = mpf_get_d (re_b);
		double ia = mpf_get_d (im_b);
		if (rt*rt +it*it < 1.0e-15 * (ra*ra+ia*ia)) 
		{
			if (downer > 5) break;
			downer ++;
		}
#endif

	}

	/* add const terms */
	mpf_set_d (term, re_s);
	mpf_mul (term, term, gam);
	mpf_sub (re_b, re_b, term);

	mpf_set_d (term, im_s);
	mpf_mul (term, term, gam);
	mpf_sub (im_b, im_b, term);

	/* subtract 1/2 */
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

/* ======================================================================= */
/* compute entire_sub_s for complex-valued s
 */

void entire_sub_s (mpf_t re_b, mpf_t im_b, double re_s, double im_s, unsigned int prec, int nterms)
{
	int k;
	mpf_t rebin, imbin, term, racc, iacc, rzeta, izeta;

	mpf_init (term);
	mpf_init (racc);
	mpf_init (iacc);
	mpf_init (rzeta);
	mpf_init (izeta);
	mpf_init (rebin);
	mpf_init (imbin);
	
	mpf_set_ui (re_b, 0);
	mpf_set_ui (im_b, 0);

	int n = 650;  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	n = nterms;
	int downer  = 0;
	for (k=2; k<= n; k++)
	{
		/* Commpute the binomial */
		c_binomial (rebin, imbin, re_s, im_s, k);

// printf ("duude s= (%g %g) k=%d bin=(%g %g)\n", re_s, im_s, k, mpf_get_d(rebin), mpf_get_d(imbin));

		/* compute zeta (k) */
		fp_zeta (term, k, prec);
		mpf_sub_ui (term, term, 1);

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
#if 0
		double rt = mpf_get_d (rzeta);
		double it = mpf_get_d (izeta);
		double ra = mpf_get_d (re_b);
		double ia = mpf_get_d (im_b);
		if (rt*rt +it*it < 1.0e-15 * (ra*ra+ia*ia)) 
		{
			if (downer > 5) break;
			downer ++;
		}
#endif

	}

	mpf_clear (term);
	mpf_clear (racc);
	mpf_clear (iacc);
	mpf_clear (rzeta);
	mpf_clear (izeta);
	mpf_clear (rebin);
	mpf_clear (imbin);
}

/* ============================================================================= */

static mpf_t re_a, im_a;
static int prec;
static int nterms;

static void a_s_init (void)
{
	/* the decimal precison (number of decimal places) */
	prec = 300;
   nterms = 300;

	/* compute number of binary bits this corresponds to. */
	double v = ((double) prec) *log(10.0) / log(2.0);

	/* the variable-precision calculations are touchy about this */
	/* XXX this should be stirling's approx for binomial */ 
	int bits = (int) (v + 30);
	
	/* Set the precision (number of binary bits) */
	mpf_set_default_prec (bits);
	
	mpf_init (re_a);
	mpf_init (im_a);
}
	
static double a_s (double re_s, double im_s)
{
	// a_sub_s (re_a, im_a, re_s, im_s, prec);
	// b_sub_s (re_a, im_a, re_s, im_s, prec, nterms);
	entire_sub_s (re_a, im_a, re_s, im_s, prec, nterms);

	double frea = mpf_get_d (re_a);
	double fima = mpf_get_d (im_a);

	double phase = atan2 (fima, frea);
	phase += M_PI;
	phase /= 2.0*M_PI;
	return phase;
}

/*-------------------------------------------------------------------*/
/* This routine fills in the interior of the the convergent area of the 
 * Euler totient in a simple way 
 */

void 
MakeHisto (
   float  	*glob,
   int 		sizex,
   int 		sizey,
   double	re_center,
   double	im_center,
   double	width,
   double	height,
   int		itermax,
   double 	renorm)
{
   int		i,j, globlen;
   double	re_start, im_start, delta;
   double	re_position, im_position;
   
   delta = width / (double) sizex;
   re_start = re_center - width / 2.0;
   im_start = im_center + width * ((double) sizey) / (2.0 * (double) sizex);
   
   globlen = sizex*sizey;
   for (i=0; i<globlen; i++) glob [i] = 0.0;

	prec = itermax;
	a_s_init();

   im_position = im_start;
   for (i=0; i<sizey; i++) 
	{
      // if (i%10==0) printf(" start row %d\n", i);
      printf(" start row %d\n", i);
      re_position = re_start;
      for (j=0; j<sizex; j++) 
		{

			double phi = a_s (re_position, im_position);
         glob [i*sizex +j] = phi;

         re_position += delta;
      }
      im_position -= delta;  /*top to bottom, not bottom to top */
   }
}

/* --------------------------- END OF LIFE ------------------------- */

/*
 * mp-cache.c
 *
 * Array cachine functions for arrays of high precision quantities
 * Gnu Multiple-precision library.
 *
 * Linas Vepstas July 2005
 */

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>
#include "db-cache.h"

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
static int i_one_d_cache_check (i_cache *c, unsigned int n)
{
	if (c->disabled) return 0;
	if ((n > c->nmax) || 0==n )
	{
		unsigned int newsize = 1.5*n+1;
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
/** i_triangle_cache_check() -- check if mpz_t value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a trianglular cache layout (two indecies)
 *  with 0 <= k <=n
 */
int i_triangle_cache_check (i_cache *c, unsigned int n, unsigned int k)
{
	if ((n > c->nmax) || 0==n )
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
	if ((n > c->nmax) || 0==n )
	{
		unsigned int newsize = 1.5*n+1;
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
	static fp_cache name = {.nmax=0, .cache=NULL, .precision=NULL}

/** fp_one_d_cache_check() -- check if mpf_t value is in the cache
 *  If there is a cached value, this returns the precision of the 
 *  value in the cache; else it returns zero.
 *  This assumes a 1-dimensional cache layout (simple array)
 */
int fp_one_d_cache_check (fp_cache *c, unsigned int n)
{
	if ((n > c->nmax) || 0==n )
	{
		unsigned int newsize = 1.5*n+1;
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
/** fp_triangle_cache_check() -- check if mpf_t value is in the cache
 *  If there is a cached value, this returns the precision of the 
 *  value in the cache; else it returns zero.
 *  This assumes a trianglular cache layout (two indecies)
 *  with 0 <= k <=n
 */
int fp_triangle_cache_check (fp_cache *c, unsigned int n, unsigned int k)
{
	if ((n > c->nmax) || 0==n )
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

/* =============================== END OF FILE =========================== */


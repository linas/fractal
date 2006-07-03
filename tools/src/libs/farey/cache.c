
/* cache.c
 * Generic cache management for commonly computed numbers
 *
 * Linas Vepstas 2005,2006
 */

#include <stdlib.h>

#include "cache.h"

/** ld_one_d_cache_check() -- check if long double value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int ld_one_d_cache_check (ld_cache *c, unsigned int n)
{
	if (c->disabled) return 0;
	if ((n > c->nmax) || 0==n )
	{
		unsigned int newsize = 1.5*n+1;
		c->cache = (long double *) realloc (c->cache, newsize * sizeof (long double));
		c->ticky = (char *) realloc (c->ticky, newsize * sizeof (char));

		unsigned int en;
		unsigned int nstart = c->nmax+1;
		if (0 == c->nmax) nstart = 0;
		for (en=nstart; en <newsize; en++)
		{
			c->cache[en] = 0.0L;
			c->ticky[en] = 0;
		}
		c->nmax = newsize-1;
		return 0;
	}

	return (c->ticky[n]);
}

/** 
 * ld_one_d_cache_fetch - fetch value from cache
 */
long double ld_one_d_cache_fetch (ld_cache *c, unsigned int n)
{
	if (c->disabled) return 0.0L;
	return c->cache[n];
}

/**
 * ld_one_d_cache_store - store value in cache
 */
void ld_one_d_cache_store (ld_cache *c, long double val, unsigned int n)
{
	if (c->disabled) return;
	c->cache[n] = val;
	c->ticky[n] = 1;
}


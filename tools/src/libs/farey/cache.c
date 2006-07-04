
/* cache.c
 * Generic cache management for commonly computed numbers
 *
 * Linas Vepstas 2005,2006
 */

#include <stdlib.h>

#include "cache.h"

/* =============================================================== */
/** TYPE_NAME##_d_cache_check() -- check if long double value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
#define CACHE_CHECK(TYPE_NAME,TYPE) \
int TYPE_NAME##_one_d_cache_check (TYPE_NAME##_cache *c, unsigned int n)	\
{	\
	if (c->disabled) return 0;	\
	if ((n > c->nmax) || 0==n )	\
	{	\
		unsigned int newsize = 1.5*n+1;	\
		c->cache = (TYPE *) realloc (c->cache, newsize * sizeof (TYPE));	\
		c->ticky = (char *) realloc (c->ticky, newsize * sizeof (char));	\
	\
		unsigned int en;	\
		unsigned int nstart = c->nmax+1;	\
		if (0 == c->nmax) nstart = 0;	\
		for (en=nstart; en <newsize; en++)	\
		{	\
			c->cache[en] = 0.0L;	\
			c->ticky[en] = 0;	\
		}	\
		c->nmax = newsize-1;	\
		return 0;	\
	}	\
	\
	return (c->ticky[n]);	\
}

/** 
 * TYPE_NAME##_d_cache_fetch - fetch value from cache	
 */
#define CACHE_FETCH(TYPE_NAME,TYPE) \
TYPE TYPE_NAME##_one_d_cache_fetch (TYPE_NAME##_cache *c, unsigned int n)	\
{	\
	if (c->disabled) return 0.0L;	\
	return c->cache[n];	\
}

/**
 * TYPE_NAME##_d_cache_store - store value in cache
 */
#define CACHE_STORE(TYPE_NAME,TYPE) \
void TYPE_NAME##_one_d_cache_store (TYPE_NAME##_cache *c, TYPE val, unsigned int n)	\
{	\
	if (c->disabled) return;	\
	c->cache[n] = val;	\
	c->ticky[n] = 1;	\
}

/* =============================================================== */

#define DEFINE_CACHE(TYPE_NAME,TYPE) \
		  CACHE_CHECK(TYPE_NAME,TYPE) \
		  CACHE_FETCH(TYPE_NAME,TYPE) \
		  CACHE_STORE(TYPE_NAME,TYPE)

DEFINE_CACHE(ld, long double)
DEFINE_CACHE(ui, unsigned int)
		  

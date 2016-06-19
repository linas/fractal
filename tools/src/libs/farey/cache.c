/**
 * cache.c
 * Generic cache management for frequently requested numbers.
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
#define CACHE_CHECK(TYPE_NAME,TYPE,IDX_TYPE) \
bool TYPE_NAME##_one_d_cache_check(TYPE_NAME##_cache *c, IDX_TYPE n)	\
{	\
	if (c->disabled) return false;	\
	if ((n > c->nmax) || 0==n )	\
	{	\
		IDX_TYPE newsize = 1.5*n+100;	\
		c->cache = (TYPE *) realloc (c->cache, newsize * sizeof (TYPE));	\
		c->ticky = (bool *) realloc (c->ticky, newsize * sizeof (bool));	\
		\
		IDX_TYPE en;	\
		IDX_TYPE nstart = c->nmax+1;	\
		if (0 == c->nmax) nstart = 0;	\
		for (en=nstart; en <newsize; en++)	\
		{	\
			c->cache[en] = (TYPE) 0;	\
			c->ticky[en] = false;	\
		}	\
		c->nmax = newsize-1;	\
		return false;	\
	}	\
	\
	return c->ticky[n];	\
}

/**
 * TYPE_NAME##_d_cache_fetch - fetch value from cache
 */
#define CACHE_FETCH(TYPE_NAME,TYPE,IDX_TYPE) \
TYPE TYPE_NAME##_one_d_cache_fetch(TYPE_NAME##_cache *c, IDX_TYPE n)	\
{	\
	if (c->disabled) return (TYPE) 0;	\
	return c->cache[n];	\
}

/**
 * TYPE_NAME##_d_cache_store - store value in cache
 */
#define CACHE_STORE(TYPE_NAME,TYPE,IDX_TYPE) \
void TYPE_NAME##_one_d_cache_store(TYPE_NAME##_cache *c, TYPE val, IDX_TYPE n)	\
{	\
	if (c->disabled) return;	\
	c->cache[n] = val;	\
	c->ticky[n] = true;	\
}

/**
 * TYPE_NAME##_d_cache_clear - clear the cache
 */
#define CACHE_CLEAR(TYPE_NAME,TYPE,IDX_TYPE) \
void TYPE_NAME##_one_d_cache_clear(TYPE_NAME##_cache *c)	\
{	\
	IDX_TYPE en; \
	for (en=0; en < c->nmax; en++)	\
	{	\
		c->cache[en] = (TYPE) 0;	\
		c->ticky[en] = false;	\
	}	\
}

/* =============================================================== */

#define DEFINE_CACHE(TYPE_NAME,TYPE,IDX_TYPE) \
		  CACHE_CHECK(TYPE_NAME,TYPE,IDX_TYPE) \
		  CACHE_FETCH(TYPE_NAME,TYPE,IDX_TYPE) \
		  CACHE_STORE(TYPE_NAME,TYPE,IDX_TYPE) \
		  CACHE_CLEAR(TYPE_NAME,TYPE,IDX_TYPE)

DEFINE_CACHE(ld, long double, unsigned int)
DEFINE_CACHE(ui, unsigned int, unsigned int)
DEFINE_CACHE(ul, unsigned long, unsigned long)

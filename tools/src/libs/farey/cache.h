/**
 * cache.h
 * Generic cache management for frequently requested numbers.
 *
 * Linas Vepstas 2005,2006
 */

#include <stdbool.h>

/* ======================================================================= */
/* Cache management */

typedef struct {
	unsigned int nmax;
	long double *cache;
	bool *ticky;
	bool disabled;
} ld_cache;


#define DECLARE_LD_CACHE(name)         \
	static ld_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, .disabled = false}

/** ld_one_d_cache_check() -- check if long double value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple array)
 */
bool ld_one_d_cache_check (ld_cache *c, unsigned int n);

/** 
 * ld_one_d_cache_fetch - fetch value from cache
 */
long double ld_one_d_cache_fetch (ld_cache *c, unsigned int n);

/**
 * ld_one_d_cache_store - store value in cache
 */
void ld_one_d_cache_store (ld_cache *c, long double val, unsigned int n);

/* ======================================================================= */
/* Cache management */

typedef struct {
	unsigned int nmax;
	unsigned int *cache;
	bool *ticky;
	bool disabled;
} ui_cache;

#define DECLARE_UI_CACHE(name)         \
	static ui_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, .disabled = false}

/** ui_one_d_cache_check() -- check if long double value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple array)
 */
bool ui_one_d_cache_check (ui_cache *c, unsigned int n);

/** 
 * ui_one_d_cache_fetch - fetch value from cache
 */
unsigned int ui_one_d_cache_fetch (ui_cache *c, unsigned int n);

/**
 * ui_one_d_cache_store - store value in cache
 */
void ui_one_d_cache_store (ui_cache *c, unsigned int val, unsigned int n);


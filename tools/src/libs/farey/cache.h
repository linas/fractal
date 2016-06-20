/**
 * cache.h
 * Generic cache management for frequently requested numbers.
 *
 * Linas Vepstas 2005,2006
 */

#include <pthread.h>
#include <stdbool.h> // for boolean
#include <stddef.h>  // for NULL

/* ======================================================================= */
/* Cache management - long double. */

typedef struct {
	unsigned int nmax;
	long double *cache;
	bool *ticky;
	bool disabled;
	pthread_rwlock_t lck;
} ld_cache;


#define DECLARE_LD_CACHE(name)         \
	static ld_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, \
		.disabled = false, .lck = PTHREAD_RWLOCK_INITIALIZER}

/** ld_one_d_cache_check() -- check if long double value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple array).
 */
bool ld_one_d_cache_check(ld_cache *c, unsigned int n);

/** 
 * ld_one_d_cache_fetch - fetch value from cache
 */
long double ld_one_d_cache_fetch(ld_cache *c, unsigned int n);

/**
 * ld_one_d_cache_store - store value in cache
 */
void ld_one_d_cache_store(ld_cache *c, long double val, unsigned int n);

/**
 * Clear the cache.
 */
void ld_one_d_cache_clear(ld_cache *c);

/* ======================================================================= */
/* Cache management -- unsigned int */

typedef struct {
	unsigned int nmax;
	unsigned int *cache;
	bool *ticky;
	bool disabled;
	pthread_rwlock_t lck;
} ui_cache;

#define DECLARE_UI_CACHE(name)         \
	static ui_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, \
		.disabled = false, .lck = PTHREAD_RWLOCK_INITIALIZER}

/** ui_one_d_cache_check() -- check if the uint value is in the cache.
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple array).
 */
bool ui_one_d_cache_check(ui_cache *c, unsigned int n);

/** 
 * ui_one_d_cache_fetch - fetch value from cache
 */
unsigned int ui_one_d_cache_fetch(ui_cache *c, unsigned int n);

/**
 * ui_one_d_cache_store - store value in cache
 */
void ui_one_d_cache_store(ui_cache *c, unsigned int val, unsigned int n);

/**
 * Clear the cache.
 */
void ui_one_d_cache_clear(ui_cache *c);

/* ======================================================================= */
/* Cache management -- unsigned long */

typedef struct {
	unsigned long nmax;
	unsigned long *cache;
	bool *ticky;
	bool disabled;
	pthread_rwlock_t lck;
} ul_cache;

#define DECLARE_UL_CACHE(name)         \
	static ul_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, \
		.disabled = false, .lck = PTHREAD_RWLOCK_INITIALIZER}

/** ul_one_d_cache_check() -- check if the ulong value is in the cache.
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple array).
 */
bool ul_one_d_cache_check(ul_cache *c, unsigned long n);

/** 
 * ul_one_d_cache_fetch - fetch value from cache
 */
unsigned long ul_one_d_cache_fetch(ul_cache *c, unsigned long n);

/**
 * ul_one_d_cache_store - store value in cache
 */
void ul_one_d_cache_store(ul_cache *c, unsigned long val, unsigned long n);

/**
 * Clear the cache.
 */
void ul_one_d_cache_clear(ul_cache *c);

/* ======================================================================= */

/**
 * cache.h
 * Generic cache management for frequently requested numbers.
 *
 * Linas Vepstas 2005,2006
 */

#include <complex.h>
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
	pthread_spinlock_t spin;
} ld_cache;


#define DECLARE_LD_CACHE(name)         \
	static ld_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, \
		.disabled = false }; \
	__attribute__((constructor)) void name##_init(void) { \
	pthread_spin_init(&name.spin, 0); }

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
/* Cache management - double complex. */

//  C and C++ is fugnuts insane in complex support.
#define complex _Complex

typedef struct {
	unsigned int nmax;
	complex *cache;
	bool *ticky;
	bool disabled;
	pthread_spinlock_t spin;
} cpx_cache;


#define DECLARE_CPX_CACHE(name)         \
	static cpx_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, \
		.disabled = false }; \
	__attribute__((constructor)) void name##_init(void) { \
	pthread_spin_init(&name.spin, 0); }

/** cpx_one_d_cache_check() -- check if double complex value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple array).
 */
bool cpx_one_d_cache_check(cpx_cache *c, unsigned int n);

/**
 * ld_one_d_cache_fetch - fetch value from cache
 */
complex cpx_one_d_cache_fetch(cpx_cache *c, unsigned int n);

/**
 * ld_one_d_cache_store - store value in cache
 */
void cpx_one_d_cache_store(cpx_cache *c, complex val, unsigned int n);

/**
 * Clear the cache.
 */
void cpx_one_d_cache_clear(cpx_cache *c);

/* ======================================================================= */
/* Cache management -- unsigned int */

typedef struct {
	unsigned int nmax;
	unsigned int *cache;
	bool *ticky;
	bool disabled;
	pthread_spinlock_t spin;
} ui_cache;

#define DECLARE_UI_CACHE(name)         \
	static ui_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, \
		.disabled = false, }; \
	__attribute__((constructor)) void name##_init(void) { \
	pthread_spin_init(&name.spin, 0); }

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
	pthread_spinlock_t spin;
} ul_cache;

#define DECLARE_UL_CACHE(name)         \
	static ul_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, \
		.disabled = false, }; \
	__attribute__((constructor)) void name##_init(void) { \
	pthread_spin_init(&name.spin, 0); }

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
/* Cache management -- unsigned long long */

typedef struct {
	unsigned __int128 nmax;
	unsigned __int128 *cache;
	bool *ticky;
	bool disabled;
	pthread_spinlock_t spin;
} ull_cache;

#define DECLARE_ULL_CACHE(name)         \
	static ull_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, \
		.disabled = false, }; \
	__attribute__((constructor)) void name##_init(void) { \
	pthread_spin_init(&name.spin, 0); }

/** ull_one_d_cache_check() -- check if the u__int128 value is in the cache.
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple array).
 */
bool ull_one_d_cache_check(ull_cache *c, unsigned int n);

/**
 * ull_one_d_cache_fetch - fetch value from cache
 */
unsigned __int128 ull_one_d_cache_fetch(ull_cache *c, unsigned int n);

/**
 * ull_one_d_cache_store - store value in cache
 */
void ull_one_d_cache_store(ull_cache *c, unsigned __int128 val, unsigned int n);

/**
 * Clear the cache.
 */
void ull_one_d_cache_clear(ull_cache *c);

/* ======================================================================= */

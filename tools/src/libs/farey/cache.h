
/* cache.h
 * Generic cache management for commonly computed numbers
 *
 * Linas Vepstas 2005,2006
 */

/* ======================================================================= */
/* Cache management */

typedef struct {
	unsigned int nmax;
	long double *cache;
	char *ticky;
	short disabled;
} ld_cache;


#define DECLARE_LD_CACHE(name)         \
	static ld_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, .disabled = 0}

/** ld_one_d_cache_check() -- check if long double value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int ld_one_d_cache_check (ld_cache *c, unsigned int n);

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
	char *ticky;
	short disabled;
} ui_cache;

#define DECLARE_UI_CACHE(name)         \
	static ui_cache name = {.nmax=0, .cache=NULL, .ticky=NULL, .disabled = 0}

/** ui_one_d_cache_check() -- check if long double value is in the cache
 *  Returns true if the value is in the cache, else returns false.
 *  This assumes a 1-dimensional cache layout (simple aray)
 */
int ui_one_d_cache_check (ui_cache *c, unsigned int n);

/** 
 * ui_one_d_cache_fetch - fetch value from cache
 */
unsigned int ui_one_d_cache_fetch (ui_cache *c, unsigned int n);

/**
 * ui_one_d_cache_store - store value in cache
 */
void ui_one_d_cache_store (ui_cache *c, unsigned int val, unsigned int n);


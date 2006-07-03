
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


/*
 * FUNCTION:
 * return greatest prime factor (larged prime divisor)
 *
 * HISTORY:
 * April 2016 -- linas
 */

#include "cache.h"
#include "gpf.h"
#include "prime.h"

/* ------------------------------------------------------------ */
/**
 * Return the greatest prime factor.
 * Direct computation (from scratch) each time -- no caching.
 */
static unsigned long gpf_direct(unsigned long n)
{
	unsigned long fact = 1;
	for (unsigned int nth = 1; ; nth++)
	{
		unsigned long p = get_nth_prime(nth);
		if (n < p) return fact;

		if (n % p == 0)
		{
			fact = p;
			n /= p;
			while (n % p == 0) n /= p;
		}
	}

	return 0;
}

/**
 * Return the greatest prime factor.
 * Cached version -- avoids recomputation.
 */
unsigned long gpf(unsigned long n)
{
	DECLARE_UL_CACHE(gpf_cache);
	if (ul_one_d_cache_check(&gpf_cache, n))
		return ul_one_d_cache_fetch(&gpf_cache, n);
	unsigned long fact = gpf_direct(n);
	ul_one_d_cache_store(&gpf_cache, fact, n);
	return fact;
}

// #define TEST 1
#ifdef TEST
#include <stdio.h>

int main()
{
	for (unsigned long n=1; n<100; n++)
	{
		printf("n=%lu gpf=%lu\n", n, gpf(n));
	}
}
#endif

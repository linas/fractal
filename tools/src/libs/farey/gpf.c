/*
 * FUNCTION:
 * return greatest prime factor (largest prime divisor)
 *
 * HISTORY:
 * April 2016 -- linas
 */

#include <pthread.h>
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

static pthread_mutex_t lck;
static void gpf_init(void)  __attribute__((constructor));
static void gpf_init(void) 
{
	pthread_mutexattr_t attr;
	pthread_mutexattr_init(&attr);
	pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_NORMAL);
	pthread_mutex_init(&lck, &attr);
	pthread_mutexattr_destroy(&attr);
}

/**
 * Return the greatest prime factor.
 * Cached version -- avoids recomputation.
 */
unsigned long gpf(unsigned long n)
{
	DECLARE_UL_CACHE(gpf_cache);
	pthread_mutex_lock(&lck);
	if (ul_one_d_cache_check(&gpf_cache, n))
	{
		unsigned long fact = ul_one_d_cache_fetch(&gpf_cache, n);
		pthread_mutex_unlock(&lck);
		return fact;
	}
	unsigned long fact = gpf_direct(n);
	ul_one_d_cache_store(&gpf_cache, fact, n);
	pthread_mutex_unlock(&lck);
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

/*
 * FUNCTION:
 * return greatest prime factor (largest prime divisor)
 *
 * HISTORY:
 * April 2016 -- linas
 */

#include <stdlib.h>

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

DECLARE_UL_CACHE(gpf_cache);

/**
 * Return the greatest prime factor.
 * Cached version -- avoids recomputation.
 */
unsigned long gpf(unsigned long n)
{
	if (ul_one_d_cache_check(&gpf_cache, n))
	{
		return ul_one_d_cache_fetch(&gpf_cache, n);
	}
	unsigned long fact = gpf_direct(n);
	ul_one_d_cache_store(&gpf_cache, fact, n);
	return fact;
}

/* ------------------------------------------------------------ */
/**
 * Return a random number with a distribution similar to the greatest
 * prime factor.
 * Direct computation (from scratch) each time -- no caching.
 */
static unsigned long pseudo_gpf_direct(unsigned long n)
{
	/* Compute the approximate distribution. It is given by the
	 * harmonic nomber of the number of primes less than n.
	 */
	double scale = 0.0;
	unsigned int nth;
	for (nth = 1; ; nth++)
	{
		unsigned long p = get_nth_prime(nth);
		if (n < p) break;
		scale += 1.0 / nth;
	}
	unsigned int num_primes = nth-1;

	/* Generate a random prime number... */
	/* Use the sliding scale from above to do so. */
	double ran = random();
	ran /= RAND_MAX;
	ran *= scale;
	double cut = 0.0;
	for (nth = 1; nth < num_primes; nth++)
	{
		cut += 1.0/nth;
		if (ran < cut) return get_nth_prime(1+num_primes-nth);
	}

	return 2; /* unlikely but can happen */
}

DECLARE_UL_CACHE(pseudo_gpf_cache);

/**
 * Return the greatest prime factor.
 * Cached version -- avoids recomputation.
 */
unsigned long pseudo_gpf(unsigned long n)
{
	if (ul_one_d_cache_check(&pseudo_gpf_cache, n))
	{
		return ul_one_d_cache_fetch(&pseudo_gpf_cache, n);
	}
	unsigned long fact = pseudo_gpf_direct(n);
	ul_one_d_cache_store(&pseudo_gpf_cache, fact, n);
	return fact;
}

/* ------------------------------------------------------------ */
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

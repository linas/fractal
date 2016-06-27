/*
 * FUNCTION:
 * return greatest prime factor (largest prime divisor)
 *
 * HISTORY:
 * April 2016 -- linas
 */

#include <math.h>
#include <stdio.h>
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
	if (n <= 1) return 1;

	/* How many primes are there less than or equal to n? */
	unsigned int nth;
	for (nth = 1; ; nth++)
	{
		unsigned long p = get_nth_prime(nth);
		if (n < p) break;
	}
	unsigned int num_primes = nth-1;

	/* Compute the approximate distribution. It is given by the
	 * tweaked harmonic number of ... plausible divisors of n...
	 */
	double punt = get_nth_prime(num_primes) / 2.0;
	double scale = 1.0;
	int cnt = 1;
	for (nth = num_primes-1; nth>0; nth--)
	{
		if (get_nth_prime(nth) > punt) continue;
		cnt ++;
		punt *= cnt;
		punt /= (cnt+1);
		double x = 1.0 / cnt;
		// scale += x;
		scale += x - 0.1*x*log(x);
	}

	/* Generate a random prime number... */
	/* Use the sliding scale from above to do so. */
	double ran = random();
	ran /= RAND_MAX;
	ran *= scale;

	double cut = 0.0;
	punt = get_nth_prime(num_primes) / 2.0;
	cnt = 1;
	for (nth = num_primes-1; nth>0; nth--)
	{
		if (get_nth_prime(nth) > punt) continue;
		cnt ++;
		punt *= cnt;
		punt /= (cnt+1);
		double x = 1.0 / cnt;
		// cut += x;
		cut += x - 0.1*x*log(x);
		if (ran < cut) break;
	}
	if (0 == nth) nth=1;

	return get_nth_prime(nth);
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
#define TEST 1
#ifdef TEST
#include <stdio.h>

int main()
{
	int hilo = 0;
	for (unsigned long n=1; n<500; n++)
	{
		if (gpf(n) < pseudo_gpf(n)) hilo++; else hilo--;
		printf("n=%lu gpf=%lu  pseudo=%lu hilo=%d\n", n, gpf(n), pseudo_gpf(n), hilo);
	}
}
#endif

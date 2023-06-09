/*
 * FUNCTION:
 * Return greatest prime factor (largest prime divisor)
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
	if (0 == n) return 0;

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
	unsigned int num_primes = prime_count(n);

	/* Compute the approximate distribution. It is given by the
	 * tweaked harmonic number of ... plausible divisors of n...
	 * This started out trying to be a scientific, rational estimate
	 * of the probability distribution, but this proved to be too
	 * complicated to get right, and devolved into a total hack.
	 * It gives a distribution that goes too low for small n,
	 * and goes too high for large n.  Its vaguely-close-ish.
	 * Kind-of-ish.
	 */
	double punt = get_nth_prime(num_primes) / 2.0;
	double scale = 1.0;
	int cnt = 1;
	unsigned int nth;
	for (nth = num_primes-1; nth>0; nth--)
	{
		if (get_nth_prime(nth) > punt) continue;
		cnt ++;
		punt *= cnt;
		punt /= (cnt+1);
		double x = 1.0 / cnt;
		// scale += x;
		scale += x - 0.1*x*log(x);
		punt /= scale;   // total bad hackery
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
		punt /= cut; // Total bad hackery ...
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
/**
 * Return the product of all of the prime factors.
 * Each factor occurs only once: squares and higher are collapsed.
 * Direct computation (from scratch) each time -- no caching.
 */
static unsigned long factor_product_direct(unsigned long n)
{
	unsigned long prod = 1;
	for (unsigned int nth = 1; ; nth++)
	{
		unsigned long p = get_nth_prime(nth);
		if (n < p) return prod;

		if (n % p == 0)
		{
			prod *= p;
			n /= p;
			while (n % p == 0) n /= p;
		}
	}

	return 0;
}

DECLARE_UL_CACHE(factor_prod_cache);

/**
 * Return the product of all of the prime factors.
 * Cached version -- avoids recomputation.
 */
unsigned long factor_product(unsigned long n)
{
	if (ul_one_d_cache_check(&factor_prod_cache, n))
	{
		return ul_one_d_cache_fetch(&factor_prod_cache, n);
	}
	unsigned long fact = factor_product_direct(n);
	ul_one_d_cache_store(&factor_prod_cache, fact, n);
	return fact;
}

/* ------------------------------------------------------------ */
// #define TEST 1
#ifdef TEST
#include <stdio.h>

int main()
{
	int hilo = 0;
	for (unsigned long n=1; n<3500; n++)
	{
		if (gpf(n) < pseudo_gpf(n)) hilo++; else hilo--;
		printf("n=%lu gpf=%lu  pseudo=%lu hilo=%d\n", n, gpf(n), pseudo_gpf(n), hilo);
	}
}
#endif

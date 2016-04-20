/*
 * FUNCTION:
 * (Very) basic prime-number sieve.
 *
 * HISTORY:
 * Linas Vepstas January 2005
 * Updates July 2006
 * Updates November 2014
 */


#include <math.h>
#include <malloc.h>
#include <pthread.h>
#include "prime.h"

/* An unsigned int32 is sufficient for factoring 64-bit ints. */
static unsigned int *sieve = NULL;
static size_t sieve_size = 0;  /* size, in bytes. */
static size_t sieve_max = 0;   /* largest correct entry. */

#define INIT_PRIME_SIEVE(N) \
	if (!sieve || sieve[sieve_max]*sieve[sieve_max] <(N)) {\
		init_prime_sieve(N); \
	}

/* Initialize and fill in a prime-number sieve.
 * Handles primes up to 4 billion (2^32), which should be enough
 * to factor 64-bit numbers.
 *
 * Initializes size so that it contains at least np primes in it.
 */
static void
init_prime_sieve (size_t max)
{
	static pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;

	unsigned int n, j;
	unsigned int nstart;
	unsigned int pos;

	if (max < sieve_max) return;

	pthread_mutex_lock(&mtx);

	// Test again, this time under the lock.
	if (max < sieve_max)
	{
		pthread_mutex_unlock(&mtx);
		return;
	}

	sieve_size = 8192 * ((max / 8192) + 1);
	if (!sieve)
	{
		sieve = (unsigned int *) malloc(sieve_size * sizeof(unsigned int));
		sieve_max = 2;
		sieve[0] = 2;
		sieve[1] = 3;
		sieve[2] = 5;
	}
	else
	{
		sieve = (unsigned int *) realloc(sieve, sieve_size * sizeof(unsigned int));
	}

	pos = sieve_max+1;
	nstart = sieve[sieve_max] + 2;

	/* Really dumb algo, brute-force test all odd numbers against
	 * known primes */
	for (n=nstart; pos <= max; n+=2)
	{
		for (j=1; ; j++)
		{
			int p = sieve[j];
			if (0 == n%p)
			{
				break;
			}
			if (p*p > n)
			{
				/* If we got to here, n must be prime; save it, move on. */
				sieve[pos] = n;
				pos ++;
				break;
			}
		}
	}
	sieve_max = pos-1;
	pthread_mutex_unlock(&mtx);

#if 0
	for (j=0; j<pos; j++)
	{
		printf ("its %d %d\n", j, sieve[j]);
	}
#endif
}

/**
 * Return the n'th prime number.
 * 2 is the 1'th prime number,
 * 3 si the 2'th prime number, etc.
 */
unsigned int get_nth_prime(unsigned long n)
{
	if (sieve_max <= n) init_prime_sieve(n);
	return sieve[n-1];
}

/* --------------------------- END OF FILE ------------------------- */

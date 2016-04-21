/*
 * FUNCTION:
 * (Very) basic prime-number sieve.
 * Thread-safe.
 *
 * HISTORY:
 * Linas Vepstas January 2005
 * Updates July 2006
 * Updates November 2014
 * Updates April 2016
 */


#include <math.h>
#include <malloc.h>
#include <pthread.h>
#include <stdbool.h>
#include <string.h>
#include "prime.h"

/* An unsigned int32 is sufficient for factoring 64-bit ints. */
static unsigned int *sieve = NULL;
static size_t sieve_size = 0;  /* size, in bytes. */
static size_t sieve_max = 0;   /* largest correct entry. */
static pthread_mutex_t mtx = PTHREAD_MUTEX_INITIALIZER;
static pthread_spinlock_t spin;

__attribute__((constructor))
void prime_sieve_spink_ctor()
{
	pthread_spin_init(&spin, 0);
}


/* Initialize and fill in a prime-number sieve.
 * Handles primes up to 4 billion (2^32), which should be enough
 * to factor 64-bit numbers.
 *
 * Initializes size so that it contains at least `max` primes in it.
 */
static void
init_prime_sieve (size_t max)
{
	if (max < sieve_max) return;

	pthread_mutex_lock(&mtx);

	// Test again, this time under the lock.
	if (max < sieve_max)
	{
		pthread_mutex_unlock(&mtx);
		return;
	}
	unsigned int pos = sieve_max+1;

	unsigned int *tsieve = sieve;
	unsigned int old_size  = sieve_size;
	sieve_size = 8192 * ((max / 8192) + 1);
	if (old_size < sieve_size)
	{
		tsieve = (unsigned int *) malloc(sieve_size * sizeof(unsigned int));
		if (!sieve)
		{
			pos = 3;
			tsieve[0] = 2;
			tsieve[1] = 3;
			tsieve[2] = 5;
		}
		else
		{
			memcpy(tsieve, sieve, old_size*sizeof(unsigned int));
		}
	}

	unsigned int nstart = tsieve[pos-1] + 2;

	/* Really dumb algo, brute-force test all odd numbers against
	 * known primes */
	for (unsigned int n=nstart; pos <= max; n+=2)
	{
		for (unsigned int j=1; ; j++)
		{
			unsigned int p = tsieve[j];
			if (0 == n%p)
			{
				break;
			}
			if (p*p > n)
			{
				/* If we got to here, n must be prime; save it, move on. */
				tsieve[pos] = n;
				pos ++;
				break;
			}
		}
	}

	pthread_spin_lock(&spin);
	if (sieve != tsieve)
	{
		unsigned int* old_sieve = sieve;
		sieve = tsieve;
		if (old_sieve) free(old_sieve);
	}
	pthread_spin_unlock(&spin);

	// Do it again, preserve memory order.
	pthread_spin_lock(&spin);
	sieve_max = pos-1;
	pthread_spin_unlock(&spin);
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

	pthread_spin_lock(&spin);
	unsigned int p = sieve[n-1];
	pthread_spin_unlock(&spin);
	return p;
}

/* --------------------------- END OF FILE ------------------------- */

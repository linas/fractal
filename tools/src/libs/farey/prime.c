/*
 * FUNCTION:
 * Basic prime-number sieve
 *
 * HISTORY:
 * Linas Vepstas January 2005
 * Updates July 2006
 * Updates November 2014
 */


#include <math.h>
#include <malloc.h>
#include "prime.h"

/* An unsigned int32 is sufficient for factoring 64-bit ints. */
static unsigned int *sieve = NULL;
static unsigned int sieve_size = 0;
static unsigned int sieve_max = 0;

#define INIT_PRIME_SIEVE(N) \
	if (!sieve || sieve[sieve_max]*sieve[sieve_max] <(N)) {\
		init_prime_sieve(N); \
	}

/* Initialize and fill in a prime-number sieve.
 * Handles primes up to 4 billion (2^32)
 * long long int should be a 64-bit number 
 */
static void
init_prime_sieve (long long int prod)
{
	unsigned int n, j;
	unsigned int nstart;
	unsigned int pos;
	
	if (sieve)
	{
		long long int ss = sieve[sieve_max];
		if (ss*ss > prod) return;
	}
		
	unsigned int max = 1000.0+sqrt (prod);

	if (!sieve)
	{
		sieve = (unsigned int *) malloc (8192*sizeof (unsigned int));
		sieve_size = 8192;
		sieve_max = 2;
		sieve[0] = 2;
		sieve[1] = 3;
		sieve[2] = 5;
	}
	pos = sieve_max+1;
	nstart = sieve[sieve_max] + 2;

	/* dumb algo, brute-force test all odd numbers against known primes */
	for (n=nstart; n<=max; n+=2)
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
				if (pos >= sieve_size)
				{
					sieve_size += 8192;
					sieve = (unsigned int *)realloc (sieve, sieve_size * sizeof (unsigned int));
				}
				break;
			}
		}
	}
	sieve_max = pos-1;

#if 0
	for (j=0; j<pos; j++)
	{
		printf ("its %d %d\n", j, sieve[j]);
	}
#endif
}

unsigned int get_nth_prime(unsigned long n)
{
	init_prime_sieve(n);
	return sieve[n];
}

/* --------------------------- END OF FILE ------------------------- */

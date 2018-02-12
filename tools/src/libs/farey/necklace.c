/*
 * necklace.c
 *
 * Return values for Moreau's necklace-counting function
 *
 * Linas Vepstas February 2018
 */

#include "moebius.h"
#include "necklace.h"
#include "cache.h"

/* Moreaus necklace counting function. OEIS A001037 */
long necklace_raw(long n)
{
	long sum = 0;
	int d = 1;
	while (d <= n/2)
	{
		if (0 == n%d) 
		{
			sum += (1<<d) * moebius_mu(n/d);
		}
		d += 1;
	}
	sum += (1<<n);

	return sum /n;
}

DECLARE_UL_CACHE (necklace_cache);

long necklace(long n)
{
	if (ul_one_d_cache_check (&necklace_cache, n))
		return ul_one_d_cache_fetch(&necklace_cache, n);

	long val = necklace_raw(n);
	ul_one_d_cache_store (&necklace_cache, val, n);
	return val;
}

// #define UNIT_TEST
#ifdef UNIT_TEST

#include <stdio.h>
int main (int arc, char* argv[])
{
	for (int i=1; i< 29; i++)
	{
		printf("%d	%ld\n", i, necklace(i));
	}
}
#endif

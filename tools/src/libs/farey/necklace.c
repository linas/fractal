/*
 * necklace.c
 *
 * Return values for Moreau's necklace-counting function
 *
 * Linas Vepstas February 2018
 */

#include <stdio.h>
#include <stdlib.h>

#include "moebius.h"
#include "necklace.h"
#include "cache.h"

/* Moreaus necklace counting function. OEIS A001037 */
long necklace_raw(int n)
{
	if (62 < n)
	{
		fprintf(stderr, "Error: necklace overflow\n");
		exit(1);
	}

	long sum = 0;
	int d = 1;
	while (d <= n/2)
	{
		if (0 == n%d) 
		{
			sum += (1UL<<d) * moebius_mu(n/d);
		}
		d += 1UL;
	}
	sum += (1UL<<n);

	return sum / n;
}

DECLARE_UL_CACHE (necklace_cache);

long necklace(int n)
{
	if (ul_one_d_cache_check (&necklace_cache, n))
		return ul_one_d_cache_fetch(&necklace_cache, n);

	long val = necklace_raw(n);
	ul_one_d_cache_store (&necklace_cache, val, n);
	return val;
}

// #define UNIT_TEST
#ifdef UNIT_TEST

int main (int arc, char* argv[])
{
	for (int i=1; i< 62; i++)
	{
		printf("%d	%ld\n", i, necklace(i));
	}
}
#endif

/*
 * dirichlet.c
 *
 * Dirichlet convolution powers of the unit function,
 * and it's inverse (the Moebius functions).
 *
 * So:
 * k=0  unit(0, 1)=1 unit(0, n)=0 for n>1
 * k=1  unit(1, n)=1 
 * k in general, unit(k,n) = (unit(k-1) * unit(1))(n)
 * k=-1 unit (-1,n) = mobius mu(n)
 * etc.
 *
 * January 2019
 */

#include <stdio.h>
#include <stdlib.h>
#include "cache.h"
#include "dirichlet.h"
#include "moebius.h"

static long unit_direct(long k, long n)
{
#if CORRECT_BUT_NOT_NEEDED
	if (n < 1) { return 0; } // actually an error.
	if (0 == k) { return n==1; }  // the "identity"
	if (1 == k) { return 1; }     // plus-one, the "unit"
	if (-1 == k) { return moebius_mu(n); } // minus-one, the inverse.
#endif // CORRECT_BUT_NOT_NEEDED

	if (1 < k)
	{
// There is a faster implementation, using the "prime cache"
// implementation used in the divisor function. But this is enough
// for just right now.
#define SLOW_BUT_SIMPLE
#ifdef SLOW_BUT_SIMPLE
		long sum = 0;
		for (long d=1; 2*d <= n; d++)
		{
			if (n%d) continue;
			sum += unit(k-1, d);
		}
		sum += unit(k-1,n);
		return sum;
#endif // SLOW_BUT_SIMPLE
	}
	if (-1 > k)
	{
#define ISLOW_BUT_SIMPLE
#ifdef ISLOW_BUT_SIMPLE
		long sum = 0;
		for (long d=1; 2*d <= n; d++)
		{
			if (n%d) continue;
			sum += unit(k+1, d) * moebius_mu(n/d);
		}
		sum += unit(k+1,n);
		return sum;
#endif // SLOW_BUT_SIMPLE
	}

	fprintf(stderr, "Internal error Dirichelt!\n");
	exit(1);
	return 0;
}

// Manual cache
// This is super-sloppy, we should use a 2D cache, but I am
// too lazy to deal with that just right now. This is enough.
DECLARE_UL_CACHE(k2_cache);
DECLARE_UL_CACHE(k3_cache);
DECLARE_UL_CACHE(k4_cache);
DECLARE_UL_CACHE(k5_cache);
DECLARE_UL_CACHE(k6_cache);
DECLARE_UL_CACHE(k7_cache);
DECLARE_UL_CACHE(k8_cache);
DECLARE_UL_CACHE(k9_cache);
DECLARE_UL_CACHE(m2_cache);
DECLARE_UL_CACHE(m3_cache);
DECLARE_UL_CACHE(m4_cache);
DECLARE_UL_CACHE(m5_cache);
DECLARE_UL_CACHE(m6_cache);
DECLARE_UL_CACHE(m7_cache);
DECLARE_UL_CACHE(m8_cache);
DECLARE_UL_CACHE(m9_cache);

#define RUN_CACHE(KAY,CAK) \
	if (KAY == k) { \
		if (ul_one_d_cache_check(&CAK, n)) \
			return ul_one_d_cache_fetch(&CAK, n); \
		long un = unit_direct(k,n); \
		ul_one_d_cache_store(&CAK, un, n); \
		return un; \
	}

long unit(long k, long n)
{
	if (n < 1) { return 0; } // actually an error.
	if (0 == k) { return n==1; }  // the "identity"
	if (1 == k) { return 1; }     // plus-one, the "unit"
	if (-1 == k) { return moebius_mu(n); } // minus-one, the inverse.

	// El slobola cacho
	RUN_CACHE(2, k2_cache);
	RUN_CACHE(-2, m2_cache);
	RUN_CACHE(3, k3_cache);
	RUN_CACHE(-3, m3_cache);
	RUN_CACHE(4, k4_cache);
	RUN_CACHE(-4, m4_cache);
	RUN_CACHE(5, k5_cache);
	RUN_CACHE(-5, m5_cache);
	RUN_CACHE(6, k6_cache);
	RUN_CACHE(-6, m6_cache);
	RUN_CACHE(7, k7_cache);
	RUN_CACHE(-7, m7_cache);
	RUN_CACHE(8, k8_cache);
	RUN_CACHE(-8, m8_cache);
	RUN_CACHE(9, k9_cache);
	RUN_CACHE(-9, m9_cache);

	return unit_direct(k,n);

	fprintf(stderr, "Internal error Dirichelt!\n");
	exit(1);
	return 0;
}

// #define TEST 1
#ifdef TEST

#include <stdio.h>

long convoid(long k, long j, long n)
{
	long sum = 0;
	for (long d=1; 2*d <= n; d++)
	{
		if (n%d) continue;
		sum += unit(k, d) * unit(j, n/d);
	}
	sum += unit(k,n) * unit(j,1);
	return sum;
}

long test_unit(void)
{
	long have_error=0;
	long nmax=10000;
	for (long n=1; n<=nmax; n++)
	{
		long convo = convoid(1, -1, n);
		if ((1 < n && 0 != convo) || (0 == n && 1 != convo))
		{
			printf ("ERROR: oh no! Moebius error at n=%ld\n", n);
			have_error ++;
		}
	}
	return have_error;
}

long test_diff(int rng, int d)
{
	long have_error=0;
	long nmax=10000;
	for (int k=-rng; k<= +rng; k++)
	{
		if (rng < d-k) continue;
		if (d-k < -rng) continue;
		for (long n=1; n<=nmax; n++)
		{
			long convo = convoid(k, d-k, n);
			long diff = unit(d, n);
			if (diff != convo)
			{
				printf ("ERROR: disadditive at diff=%d %d n=%ld got=%ld expected=%ld\n",
					k, k+d, n, convo, diff);
				have_error ++;
			}
		}
	}
	return have_error;
}


int main()
{
   long nerr = test_unit();

	for (int d=-9; d<10; d++)
	{
		nerr += test_diff(9, d);
printf("done w/ d=%d\n", d);
	}

	if (0 == nerr) printf("Dirichlet unit test passed\n");
	return nerr;
}

#endif /* TEST */

/* --------------------------- END OF FILE ------------------------- */

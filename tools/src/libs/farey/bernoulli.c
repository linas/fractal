/*
 * bernoulli.c
 *
 * Miscellaneous bernoulli-number related library routines.
 *
 * Linas Vepstas <linas@linas.org> Dec 2003, Dec 2004
 */


#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "bernoulli.h"

// ======================================================
// return the nth' bernoulli number
// this is a fairly fast algorithm ...
// must have n>=0
long double bernoulli (int n)
{
	static int nmax = 0;
	static int ncache = 0;
	static long double *bern_cache = NULL; // cache of precomputed values

	if (1 > n) return 1.0;
	if (1 == n) return -0.5;
	if (n%2 == 1) return 0.0;
	
	if (n > ncache)  // do we need more memory for the cache?
	{
		ncache = 2*n+500;
		bern_cache = realloc(bern_cache, ncache*sizeof (long double));
	}

	// is the requested value higher than what we have in cache?
	if (n > nmax) 
	{
		if (0 == nmax)
		{
			bern_cache[0] = 1.0L;
			bern_cache[1] = -0.5L;
			bern_cache[2] = 1.0L/6.0L;
			nmax = 2;
		}
		int k;

		// use recursion relation on bernoulli numbers
		for (k=nmax+2; k<n+1; k+=2) 
		{
			long double acc = 0.0L;
			acc =  bern_cache[0];
			acc += ((long double) (k+1)) * bern_cache[1];
			acc += 0.5L * ((long double)(k+1)*k) * bern_cache[2];
			int j;
			for (j=4; j<k-1; j+=2)
			{
				acc += binomial (k+1, j) * bern_cache[j];
			}
			acc /= (long double) (k+1);
			bern_cache[k] = -acc;
		}

		nmax = n; // note nmax always even here
	}

	return bern_cache[n];
}

#if TEST
int test_bernoulli (void)
{
	int i;
#define L_PI 3.1415926535897932384626433832795029L
	long double lp = 2.0L * L_PI;
	lp = logl (lp);
	long double l2 = logl (2.0L);
	int passtest = 1;
	for (i=2; i<150; i+=2)
	{
		long double bound = 2*factorial(i) / expl (((long double) i) * lp);
		long double upper = expl (((long double)(1-i)) * l2);
		upper = 1.0L/(1.0L-upper);
		upper *= bound;
		long double bern = bernoulli (i);
		if (i<50) {
			printf ("%d bernoulli=%Lf upper bound = %Lf lower bound=%Lf \n", 
							 i, bern, upper, bound);
		}
		if (i%4 == 0) bern = -bern;
		if (upper <= bern) 
		{
			printf ("ERROR: fail upper bound at i=%d\n", i);
			passtest = 0;
		}
		if (bound >= bern) 
		{
			printf ("ERROR: fail lower bound at i=%d\n", i);
			passtest = 0;
		}
	}
	if (passtest) printf ("Bernoulli test passed \n");
	printf ("\n");

	return !passtest;
}
#endif 

// ======================================================
// test harness
// 
#if TEST
main ()
{
	int failed_tests = 0;
	int total_tests = 0;
	failed_tests += test_bernoulli ();  total_tests ++;

	if (failed_tests)
	{
		printf ("Error: %d of %d tests failed\n", failed_tests, total_tests);
	}
	else
	{
		printf ("All %d tests passed\n", total_tests);
	}
}
#endif 

// ======================= END OF FILE ===================

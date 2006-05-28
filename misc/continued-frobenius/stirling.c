
/* 
 * stirling.c
 *
 * recursive stirling numbers of first and second kind 
 * 
 * Linas April 2006
 */

#include <stdio.h>
#include <stdlib.h>

#include "binomial.h"
#include "ache.h"

/** stirling_first - Return stirling numbers of the first kind.
 *
 * Return stirling numbers of the first kind,
 * normalized so that they are all positive.
 * Uses dynamically-sized cache.
 */
long double stirling_first (unsigned int n, unsigned int k)
{

	/* Cache management */
	static int nmax = 0;
	static long double *cache = NULL;
	if (n> nmax)
	{
		int newsize = n*(n+1)/2;
		cache = (long double *) realloc (cache, newsize * sizeof (long double));

		int en;
		for (en=nmax+1; en <=n; en++)
		{
			int j;
			int idx = en * (en-1) /2 - 1;
			for (j=1; j<=en; j++)
			{
				cache[idx+j] = 0.0;
			}
		}

		nmax = n;
	}

	/* Trivial case (not in the cache) */
	if (0==k)
	{
		if (0==n) return 1.0L;
		return 0.0L;
	}

	if (n<k) return 0.0L;
	if (n==k) return 1.0L;

	/* Pull value from cache if it is there */
	int idx = n * (n-1) / 2 -1;
	if (0.0 < cache[idx+k]) return cache[idx+k];

	/* use recursion to get new value */
	long double s = stirling_first (n-1, k-1);
	if (n-1 >= k)
	{
		s += (n-1) * stirling_first (n-1, k);
	}

	cache[idx+k] = s;
	return s;
}

long double gee_nk (unsigned int n, unsigned int k)
{
	/* Cache management */
	static int nmax = 0;
	static long double *cache = NULL;
	if (n> nmax)
	{
		int newsize = n*(n+1)/2;
		cache = (long double *) realloc (cache, newsize * sizeof (long double));

		int en;
		for (en=nmax+1; en <=n; en++)
		{
			int j;
			int idx = en * (en-1) /2 - 1;
			for (j=1; j<=en; j++)
			{
				cache[idx+j] = 0.0;
			}
		}

		nmax = n;
	}

	/* Trivial case (not in the cache) */
	if (0==k)
	{
		if (0==n) return -1.0L;
		if (1==n) return 1.0L;
		return 0.0L;
	}

	if (n<k) return 0.0L;
	if (n==k) return -1.0L;

	/* Pull value from cache if it is there */
	int idx = n * (n-1) / 2 -1;
	if (0.0 != cache[idx+k]) return cache[idx+k];

	/* use recursion to get new value */
	long double s = gee_nk (n-1, k-1);
	if (n-1 >= k)
	{
		s += (n-2) * gee_nk (n-1, k);
	}

	cache[idx+k] = s;
	return s;
}

long double sb_sum (unsigned int n, unsigned int m)
{
	unsigned int k;
	int sg;
	long double s;
	
	sg = 1;
	if (m%2) sg = -1;

	s = 0.0L;
	for (k=m; k<=n; k++)
	{
		long double stir = stirling_first (n,k);
		long double bin = binomial (k,m);
// printf ("duude m=%d k=%d n=%d st=%Lg  bi=%Lg\n", m, k, n, stir, bin);
		s += sg * bin * stir;
		sg = -sg;
	}

	return s;
}

long double b_sum (int n)
{
	int k;

	long double sum = 0.0;
	for (k=n; k<40; k++)
	{
		long double term = stirling_first (k,n) / factorial (k);
		term *= b_sub_n (k);
		sum += term;
	}

	return sum;
}

long double stieltjes_gamma (int n)
{
	int k;

	long double sum = 0.0L;
	long double sg = 1.0L;
	if (n%2) sg = -1.0L;
	for (k=n; k<40; k++)
	{
		long double term = b_sub_n (k);
		term *= sb_sum (k,n);
		term /= factorial (k);
		// term *= sg;
		sum += term;
		sg = -sg;
	}
	sum *= factorial (n);
	if (n%2) sum = -sum;

	return sum;
}

int
main (int argc, char *argv[]) 
{
	int n, k;

#if 1
	for (n=0; n<258; n++)
	{
		long double orm=1.0;
		int kmax = n;
		// if (10<kmax) kmax = 10;
		for (k=0; k<=kmax; k++)
		{
			// long double s = stirling_first (n,k);
			long double s = sb_sum (n,k);
			long double g = gee_nk (n,k);
			if (k%2) s=-s;
			// if (n%2) s=-s;
			if (1==k) orm = s;
			// s /= orm;
			// s /= factorial (n);
			// s *= factorial (k);
			printf ("duude (%d %d)  = %Lg  %Lg\n", n, k, s, s+g);
		}
	}
#endif

	if (2>argc)
	{
		fprintf (stderr, "Usage: %s <order>\n", argv[0]);
		exit (1);
	}

	int order = atoi (argv[1]);

	n=order;
#if 0
	for (k=0; k<=n; k++)
	{
		long double s = sb_sum (n,k);
		s /= factorial (n);
		double x = ((double) k)/ ((double) n);
		printf ("%d	%g	%Lg\n", k, x, s);
	}
#endif

#if 0
	k = order;
	for (n=k; n<k+140; n++)
	{
		long double s = sb_sum (n,k);
		// s /= factorial (n);
		// s *= factorial (k);
		printf ("%d	%Lg\n", n, s);
	}
#endif

#if 0
	for(k=0; k<40; k++)
	{
		long double val = b_sum (k);
		val *= factorial (k);
		printf ("%d	%Lg\n", k, val);
	}
#endif 

#if 0
	for (n=0; n<40; n++)
	{
		long double s = stieltjes_gamma (n);
		printf ("%d	gam= %20.16Lg\n", n, s);
	}
#endif


	return 0;
}

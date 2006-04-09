
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

long double stirling_first (unsigned int n, unsigned int k)
{
	static int nmax = 0;
	static long double *cache = NULL;
	if (n> nmax)
	{
		int newsize = n*(n+1)/2;
printf ("duude c=%p\n", cache);
		cache = (long double *) realloc (cache, newsize * sizeof (long double));

		int en;
		for (en=nmax+1; en <=n; en++)
		{
			int j;
			int idx = en * (en-1) /2 - 1;
			for (j=0; j<en; j++)
			{
				cache[idx+j] = 0.0;
			}
		}

		nmax = n;
	}

	if (0==k)
	{
		if (0==n) return 1.0L;
		return 0.0L;
	}

	int idx = n * (n-1) / 2 -1;
	if (cache[idx+k] > 0.0) return cache[idx+k];

	if (n<k) return 0.0L;

	long double s = stirling_first (n-1, k-1);
	if (n-1 >= k)
	{
		s += (n-1) * stirling_first (n-1, k);
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
	if (n%2) sg = -1;

	s = 0.0L;
	for (k=m; k<=n; k++)
	{
		s += sg * binomial (k,m) * stirling_first (n, k);
		sg = -sg;
	}

	return s;
}

int
main () 
{
	int n, k;

#if 1
	for (n=0; n<10; n++)
	{
		for (k=0; k<=n; k++)
		{
			long double s = stirling_first (n,k);
			// long double s = sb_sum (n,k);
			// s /= factorial (n);
			// s *= factorial (k);
			printf ("duude (%d %d)  = %Lg\n", n, k, s);
		}
	}
#endif

#if 0
	n=20;
	for (k=0; k<=n; k++)
	{
		long double s = sb_sum (n,k);
		s /= factorial (n);
		double x = ((double) k)/ ((double) n);
		printf ("%d	%g	%Lg\n", k, x, s);
}
#endif

	return 0;
}

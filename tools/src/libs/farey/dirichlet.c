/*
 * dirichlet.c
 *
 * Dirichlet convolution powers of the unit function,
 * and it's inverse (the Moebius functions).
 *
 * So:
 * k=0  unit(0, 1)=1 unit(0, n)=0 for n>1
 * k=1  unit(1, n)=1 
 * k in general, unit(k,n) = (unit(k) * unit(1))(n)
 * k=-1 unit (-1,n) = mobius mu(n)
 * etc.
 *
 * January 2019
 */

#include "cache.h"
#include "moebius.h"

long unit(long k, long n)
{
	if (0 == k) { return n==0; }
	if (1 == k) { return 1; }
	if (-1 == k) { return moebius_mu(n); }

	return 0;
}

#define TEST 1
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
			printf ("ERROR: oh no! n=%ld\n", n);
			have_error ++;
		}
	}
	return have_error;
}

int main()
{
   long nerr = test_unit();

	if (0 == nerr) printf("Dirichlet unit test passed\n");
	return nerr;
}

#endif /* TEST */

/* --------------------------- END OF FILE ------------------------- */

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
	if (-1 == k) {return moebius_mu(n); }
}

#define TEST 1
#ifdef TEST

#include <stdio.h>

long test_unit(void)
{
}

int main()
{
   test_unit();

	return 1;
}

#endif /* TEST */

/* --------------------------- END OF FILE ------------------------- */

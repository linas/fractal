/*
 * multi.C
 *
 * Exploration of zeta series for completely multiplicative functions
 * Linas Vepstas March 2019
 */
#include <math.h>
#include <stdio.h>

#include "euler.h"

// When this is called, p is guaranteed to be prime
static double complex at_prime (unsigned int p)
{
printf("duuude p=%d\n", p);
	double pr = (double) p;
	return sqrt(pr*pr + pr);
}

// When "arithmetic fun" is called, n is guaranteed to be prime
static double complex plicplic(arithmetic, unsigned int, unsigned int);
double complex multiplicative(arithmetic fun, unsigned int n)
{
	/* handle trivial boundary case */
	if (n <= 3) return fun(n);
	return plicplic(fun, n, 2);
}

// Helper function for above.
static complex plicplic(arithmetic fun, unsigned int y, unsigned int x)
{
	// y is prime.
	if (x+1 == y) return fun(y);

	// y is not prime
	if ((y%x) == 0)
	{
		unsigned int z = y/x;
		return multiplicative(fun, x) * multiplicative(fun, z);
	}
	else return plicplic(fun, y, x+1);
}

int main()
{
	for (int n=1; n<33; n++)
	{
		complex v = multiplicative(at_prime, n);
		printf("duuude at %d its %f+i%f\n", n, creal(v), cimag(v));
	}
}

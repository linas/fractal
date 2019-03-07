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
	double pr = (double) p;
	return sqrt(pr*pr + pr);
	// return sqrt(pr*pr - pr);
}

// Helper function
static double complex plicplic(arithmetic, unsigned int, unsigned int);

/// Return the value of a completely multiplicative function at
/// the positive integer `n`.  It is assumed that `fun` will provide
/// the values at `n`==prime. All that this function does is to
/// factor `n` into a product of primes, and then call `fun` on the
/// resulting product.
///
/// When "arithmetic fun" is called, n is guaranteed to be prime
//
double complex multiplicative(arithmetic fun, unsigned int n)
{
	/* handle trivial boundary case */
	if (n <= 3) return fun(n);
	return plicplic(fun, n, 2);
}

// Tail-recursive helper function for above.
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

void mktable()
{
	for (int n=1; n<633; n++)
	{
		complex v = multiplicative(at_prime, n);
		printf("%d	%g	%g\n", n, creal(v), cimag(v));
	}
}
int main()
{
	mktable();
}

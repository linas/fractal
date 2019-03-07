/*
 * multi.C
 *
 * Exploration of zeta series for completely multiplicative functions
 * Linas Vepstas March 2019
 */
#include <math.h>
#include <stdio.h>

#include "cache.h"
#include "euler.h"

// When this is called, p is guaranteed to be prime
static complex at_prime (unsigned int p)
{
	double pr = (double) p;
	// return sqrt(pr*pr + pr);
	// return sqrt(pr*pr - 0.5*pr);
	// return pr; // Riemann zeta
	// return sqrt(pr*pr - 0.0001*pr);
	// return sqrt(pr*pr); // Riemann again
	return sqrt(pr*pr - 1.0e-6*pr);
}

// Helper function
static complex plicplic(arithmetic, unsigned int, unsigned int);

/// Return the value of a completely multiplicative function at
/// the positive integer `n`.  It is assumed that `fun` will provide
/// the values at `n`==prime. All that this function does is to
/// factor `n` into a product of primes, and then call `fun` on the
/// resulting product.
///
/// When "arithmetic fun" is called, n is guaranteed to be prime
//
complex multiplicative(arithmetic fun, unsigned int n)
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

// Suitable for graphing with gnuplot
void mktable()
{
	for (int n=1; n<633; n++)
	{
		complex v = multiplicative(at_prime, n);
		printf("%d	%g	%g\n", n, creal(v), cimag(v));
	}
}

// Cache the values, for compute speed.
DECLARE_CPX_CACHE(plic)
complex plic_fun(unsigned int n)
{
	if (cpx_one_d_cache_check(&plic, n))
		return cpx_one_d_cache_fetch(&plic, n);

	complex val = multiplicative(at_prime, n);
	cpx_one_d_cache_store(&plic, val, n);
// printf("plic n=%d val=%f+i%f\n", n, creal(val), cimag(val));
	return val;
}

complex ess = 0.0;
complex alter_raw(unsigned int n)
{
// printf("duuude alter n=%d\n", n);
	complex val = plic_fun(n);
	val = cpow(val, -ess);
	if (n%2 == 1) return val;
	return -val;
}

DECLARE_CPX_CACHE(altern)
complex alter(unsigned int n)
{
	if (cpx_one_d_cache_check(&altern, n))
		return cpx_one_d_cache_fetch(&altern, n);

	complex val = alter_raw(n);
	cpx_one_d_cache_store(&altern, val, n);
// printf("alter n=%d val=%f+i%f\n", n, creal(val), cimag(val));
	return val;
}

int main()
{
	// mktable();
	for (double y = 0.0; y<30.0; y+=0.1)
	{
		cpx_one_d_cache_clear(&altern);
		ess = 0.5 + I*y;
		complex eta = euler_sum_cut(alter, 2500);
		complex cyc = 1.0 / (1.0 + 2.0* alter(2));
		eta *= cyc;
		printf("%g	%g	%g	%g\n", creal(ess), cimag(ess), creal(eta), cimag(eta));
		fflush(stdout);
	}
}

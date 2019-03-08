/*
 * plic.C
 *
 * Exploration of zeta series for completely multiplicative functions.
 * This is the GMP version
 *
 * Linas Vepstas March 2019
 */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mp-multiplicative.h"
#include "mp-euler.h"

complex ess = 0.0;

// When this is called, p is guaranteed to be prime
void at_prime(cpx_t f, unsigned long p, int nprec)
{
	double pr = (double) p;
	// return sqrt(pr*pr + pr);
	// return sqrt(pr*pr - 0.5*pr);
	// return sqrt(pr*pr - 0.0001*pr);
	// return sqrt(pr*pr); // Riemann again
	// return sqrt(pr*pr - 1.0e-6*pr);
	// return sqrt(pr*pr + pr);
	// return sqrt(pr*pr + 1.0e-6*pr);
	// return sqrt(pr*pr + 1.0e-4*pr);
	// return sqrt(pr*pr + 0.01*pr);

	complex val = pr; // Riemann zeta

	val = cpow(val, -ess);
	cpx_set_d(f, creal(val), cimag(val));
}

// Suitable for graphing with gnuplot
void mktable(int nprec)
{
	int nbits = 3.3*nprec;
	cpx_t result;
	cpx_init2(result, nbits);
	for (int n=1; n<633; n++)
	{
		cpx_multiplicative_cached(result, at_prime, n, nprec);
		printf("%d	%g	%g\n", n, 1.0/cpx_get_re(result), cpx_get_im(result));
	}
	cpx_clear(result);
}


#if 0
complex alter(unsigned int n)
{
// printf("duuude alter n=%d\n", n);
	complex val = plic_fun(n);
	val = cpow(val, -ess);
	if (n%2 == 1) return val;
	return -val;
}
#endif

int main()
{
	int prec, nbits;
	prec = 120;
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits+200);

#ifdef TABLE_TEST
	ess = 2.0;
	mktable(prec);
#endif

	cpx_t result;
	cpx_init2(result, nbits);

	// ess = 0.5 + I*28;
	// euler_sum_cut(alter, 2500);
exit(0);
#if 0
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
#endif
}

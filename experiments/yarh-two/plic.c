/*
 * plic.c
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
#include "mp-trig.h"

complex ess = 0.0;
cpx_t mpess;

// When this is called, p is guaranteed to be prime
void at_prime(cpx_t f, unsigned long p, int nprec)
{
	// return sqrt(pr*pr + pr);
	// return sqrt(pr*pr - 0.5*pr);
	// return sqrt(pr*pr - 0.0001*pr);
	// return sqrt(pr*pr); // Riemann again
	// return sqrt(pr*pr - 1.0e-6*pr);
	// return sqrt(pr*pr + pr);
	// return sqrt(pr*pr + 1.0e-6*pr);
	// return sqrt(pr*pr + 1.0e-4*pr);

// #define DBL_PREC
#ifdef DBL_PREC
	double pr = (double) p;
	pr = sqrt(pr*pr + 0.01*pr);
	complex val = pr; // Riemann zeta

	val = cpow(val, -ess);
	cpx_set_d(f, creal(val), cimag(val));
#endif

// #define EXACTITUDE
#ifdef EXACTITUDE
	cpx_set_ui(f, p, 0);
	cpx_recip(f, f);
	cpx_sqrt(f, f, nprec);
printf("duuude p=%lu ret=%f+i%f vs val=%f+i%f\n", p,
cpx_get_re(f), cpx_get_im(f),
creal(val), cimag(val));
#endif

	// cpx_set_ui(f, p, 0);
	double pr = (double) p;
	// pr = sqrt(pr*pr + 0.01*pr);
	// pr *= 1.01;
	// pr += 1.0;
	pr *= 1.5;
	cpx_set_d(f, pr, 0.0);
	cpx_recip(f, f);
	cpx_pow(f, f, mpess, nprec);

	// printf("duuude p=%lu ret=%f+i%f\n", p, cpx_get_re(f), cpx_get_im(f));
}

// ------------------------------------
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

// ------------------------------------
void alter(cpx_t f, unsigned long n, int nprec)
{
	cpx_multiplicative_cached(f, at_prime, n, nprec);
	if (n%2 == 1) return;
	cpx_neg(f, f); // flip the sign, alternating series
}


int main()
{
	int prec, nbits;
	prec = 320;
	nbits = 3.3*prec;
	mpf_set_default_prec (nbits+150);

	cpx_t result;
	cpx_init2(result, nbits);
	cpx_init2(mpess, nbits);

	ess = 0.5 + I*28;
	ess = 0.8 + I*30;
	ess = 0.3 + I*40;
	ess = 1.0;
	cpx_set_d(mpess, creal(ess), cimag(ess));

// #define TABLE_TEST
#ifdef TABLE_TEST
	mktable(prec);
#endif

#define SINGLE_POINT
#ifdef SINGLE_POINT
	unsigned int nterms = cpx_euler_sum(result, alter, 240, 300000, prec);

	printf("nterms=%d val=%f+i%f\n", nterms,
		cpx_get_re(result), cpx_get_im(result));
#endif
}

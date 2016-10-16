/*
 * Generating functions for assorted number-theoretic functions
 * Implementation in bignums.
 *
 * Linas Vepstas - April 2016, October 2016
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include <mp-trig.h>
#include <mp-complex.h>

#include <binomial.h>
#include <gpf.h>
#include <prime.h>
#include "gpf-gen-bignum.h"

/*
 * Ordinary generating function for the greatest prime factor.
 */
void cpx_gpf_ordinary(cpx_t sum, cpx_t z, int prec)
{
	mpf_t zabs, gabs, epsi;
	mpf_init (gabs);
	mpf_init (zabs);
	mpf_init (epsi);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.321*prec));

	cpx_set_ui(sum, 0, 0);

	// falls apart if z is zero.
	cpx_abs(gabs, z);
	if (0 > mpf_cmp(gabs, epsi)) return;

	// Not defined for |z| > 1
	mpf_sub_ui(gabs, gabs, 1);
	mpf_neg(gabs, gabs);
	if (0 > mpf_cmp(gabs, epsi)) return;

	double dist_to_circle = mpf_get_d(gabs);
	int niter = ceil (2.302585*prec / dist_to_circle);
	niter += ceil (log(niter) / dist_to_circle); // gpf bounded by n

	cpx_t zn, term;
	cpx_init(zn);
	cpx_init(term);
	cpx_set(zn, z);

	for (int n=1; n < niter ; n++)
	{
		cpx_times_ui(term, zn, gpf(n));
		cpx_add(sum, sum, term);
		cpx_mul(zn, zn, z);

#if SLOW_VERSION_NOT_USING_NITER
		// The following checks the loop termination condition,
		// which is that the size of the term is less than epsilon.
		cpx_abs(gabs, zn);
		mpf_mul_ui(gabs, gabs, n);

		cpx_abs(zabs, sum);
		mpf_div(gabs, gabs, zabs);

		// if (n * zn < epsi * sum) return;
		if (0 > mpf_cmp(gabs, epsi)) return;
#endif
	}
}


/**
 * Exponential generating function for the greatest prime factor.
 * shifted .. aka the 'offset'th derivative.
 */
void cpx_gpf_exponential_shift(cpx_t sum, cpx_t z, int offset, int prec)
{
	mpf_t zabs, gabs, epsi, fact;
	mpf_init (gabs);
	mpf_init (zabs);
	mpf_init (epsi);
	mpf_init (fact);
	mpf_set_ui(fact, 1);
	mpf_set_ui(epsi, 1);
	mpf_div_2exp(epsi, epsi, (int)(3.321*prec));

	cpx_set_ui(sum, 0, 0);

	// falls apart if z is zero.
	cpx_abs(gabs, z);
	if (0 > mpf_cmp(gabs, epsi)) return;

	cpx_t zn, term;
	cpx_init(zn);
	cpx_init(term);
	cpx_set(zn, z);

	for (int n=1; ; n++)
	{
		cpx_times_ui(term, zn, gpf(n+offset));
		cpx_times_mpf(term, term, fact);
		cpx_add(sum, sum, term);

		// The following checks the loop termination condition,
		// which is that the size of the term is less than epsilon.
		cpx_abs(gabs, term);
		mpf_mul_ui(gabs, gabs, n);

		cpx_abs(zabs, sum);
		mpf_mul(zabs, zabs, epsi);

		// if (n * zn/n! < epsi * sum) return;
		if (0 > mpf_cmp(gabs, zabs)) break;

		cpx_mul(zn, zn, z);
		mpf_div_ui(fact, fact, n+1);
	}

	// The offset is the n'th derivative.  That means that
	// there is a constant term; the above loop did not handle it.
	// add that constant now.
	if (offset != 0)
	{
		cpx_add_ui(sum, sum, gpf(offset), 0);
	}

	// Remove the leading exponential order.
	cpx_abs(gabs, z);
	mpf_neg(gabs, gabs);
	fp_exp(gabs, gabs, prec);

	cpx_times_mpf(sum, sum, gabs);
}

void cpx_gpf_exponential(cpx_t sum, cpx_t z, int prec)
{
	cpx_gpf_exponential_shift(sum, z, 0, prec);
}

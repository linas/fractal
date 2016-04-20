/*
 * Generating functions for greatest prime factors.
 * Implementation in bignums.
 *
 * April 2016
 */

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#include <mp-trig.h>
#include <mp-complex.h>

#include <gpf.h>
#include "gpf-gen-bignum.h"

/*
 * Ordinary generating function for the greatest common factor.
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
	if (0 < mpf_cmp(gabs, epsi)) return;

	cpx_t zn, term;
	cpx_init(zn);
	cpx_init(term);
	cpx_set(zn, z);

	for (int n=1; ; n++)
	{
		cpx_times_ui(term, zn, gpf(n));
		cpx_add(sum, sum, term);
		cpx_mul(zn, zn, z);

		// The following check the loop termination condition,
		// which is that the size of the term is less than epsilon.
		cpx_abs(gabs, zn);
		mpf_mul_ui(gabs, gabs, n);

		cpx_abs(zabs, sum);
		mpf_div(gabs, gabs, zabs);

		// if (n * zn < epsi * sum) return;
		if (0 > mpf_cmp(gabs, epsi)) return;
	}
}

/*
 * Exponential generating function for the greatest common factor.
 */
void cpx_gpf_exponential(cpx_t sum, cpx_t z, int prec)
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
		cpx_times_ui(term, zn, gpf(n));
		cpx_times_mpf(term, term, fact);
		cpx_add(sum, sum, term);

		// The following check the loop termination condition,
		// which is that the size of the term is less than epsilon.
		cpx_abs(gabs, term);
		mpf_mul_ui(gabs, gabs, n);

		cpx_abs(zabs, sum);
		mpf_div(gabs, gabs, zabs);

		// if (n * zn * fact < epsi * sum) return;
		if (0 > mpf_cmp(gabs, epsi)) break;

		cpx_mul(zn, zn, z);
		mpf_div_ui(fact, fact, n);
	}

	// Remove the leading exponential order.
	cpx_abs(gabs, z);
	mpf_neg(gabs, gabs);
	fp_exp(gabs, gabs, prec);

	cpx_times_mpf(sum, sum, gabs);
}

/*
 * Pochhammer generating function for the greatest common factor.
 * This is like the exponential generating function, but uses pochhammer
 * instead of x^n.  Actually, uses binomial coefficient ... two
 * factorials are needed for convergence.
 */
void cpx_gpf_poch(cpx_t sum, cpx_t z, int prec, bool rise)
{
	mpf_t zabs, gabs, fact, epsi, one;
	mpf_init (gabs);
	mpf_init (zabs);

	mpf_init (fact);
	mpf_set_ui(fact, 1);

	mpf_init (one);
	mpf_set_ui(one, 1);
	if (!rise)
		mpf_neg(one, one);

	mpf_init (epsi);
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
		cpx_times_ui(term, zn, gpf(n));
		cpx_times_mpf(term, term, fact);
		cpx_times_mpf(term, term, fact); // A second factorial!
		cpx_add(sum, sum, term);

		// The following check the loop termination condition,
		// which is that the size of the term is less than epsilon.
		cpx_abs(gabs, term);
		mpf_mul_ui(gabs, gabs, n);

		cpx_abs(zabs, sum);
		mpf_div(gabs, gabs, zabs);

		// if (n * zn/n! < epsi * sum) return;
		if (0 > mpf_cmp(gabs, epsi)) return;

		cpx_add_mpf(z, z, one);
		cpx_mul(zn, zn, z);
		mpf_div_ui(fact, fact, n);
	}
}

void cpx_gpf_poch_rising(cpx_t sum, cpx_t z, int prec)
{
	cpx_gpf_poch(sum, z, prec, true);
}

void cpx_gpf_poch_falling(cpx_t sum, cpx_t z, int prec)
{
	cpx_gpf_poch(sum, z, prec, false);
}

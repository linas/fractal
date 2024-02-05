/*
 * bigseries.c
 *
 * Compute the convergencts E_k using bignum.
 * Feb 2024
 */

#define _GNU_SOURCE
#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <anant/mp-complex.h>

/**
 * Compute binary digit sequence for a given beta value.
 *
 * digs: pointer to byte array in which to store bit sequence.
 * maxn: how many digits to compute.
 *
 * Each digit records whether 0.5<=midp (true) or not (false).
 * This is the Renyi-Parry sequence. The zeroth bit is always one.
 * the first bit is 0.5 <= (m_1 = T(beta/2))
 */
void
gen_bitseq(mpf_t beta, char* digs, int maxn)
{
	mpf_t half;
	mpf_init(half);
	mpf_set_ui(half, 1);
	mpf_div_ui(half, half, 2);

	mpf_t midp;
	mpf_init(midp);
	mpf_set(midp, beta);
	mpf_div_ui(midp, midp, 2);

	/* OK, now start iterating the map */
	for (int i=0; i < maxn; i++)
	{
		if (0 <= mpf_cmp(midp, half))
		{
			mpf_sub(midp, midp, half);
			digs[i] = 1;
		}
		else digs[i] = 0;
		mpf_mul(midp, midp, beta);
	}

	mpf_clear(half);
	mpf_clear(midp);
}

// -------------------------------------------------------
/**
 * Compute zeta polynomial. Roots of this poly will be the desired
 * zeros that we hope are accumulating.  This is the E(beta;z) thing.
 *
 * sum: returned value
 * zeta: location at which to evaluate.
 * digs: arry of digits to use
 * order: max order of the poly. There must be at least this many
 *     binary digits available. Degree of polynomial is one more than
 *     this. Sorry for the off-by-one, but so it goes.
 */
void ebz(cpx_t sum, cpx_t zeta, char* digs, int order)
{
	cpx_t zetan;
	cpx_init(zetan);
	cpx_set(zetan, zeta);

	cpx_set_ui(sum, 0, 0);

	// Do the first k-1 bits.  The last bit is always one.
	for (int i=0; i<order; i++)
	{
		if (digs[i]) cpx_add(sum, sum, zetan);
		cpx_mul(zetan, zetan, zeta);
	}
	// The final bit is always one.
	cpx_add(sum, sum, zetan);

	// subtract one.
	cpx_sub_ui(sum, sum, 1, 0);

	cpx_clear(zetan);
}

// -------------------------------------------------------
/**
 * Same as above, but computes the m'th derivative.
 * Needed as a test function for zero finding.
 */
void ebz_deriv(cpx_t sum, int mderiv, cpx_t zeta, char* digs, int order)
{
	cpx_set_ui(sum, 0, 0);
	if (order+1 < mderiv) return;

	mpf_t fact;
	mpf_init(fact);
	mpf_set_ui(fact, 1);

	cpx_t zetan;
	cpx_init(zetan);

	if (0 == mderiv)
		cpx_set(zetan, zeta);
	else
		cpx_set_ui(zetan, 1, 0);

	cpx_t term;
	cpx_init(term);

	// Leading factorial
	for (int i=1; i<mderiv; i++)
		mpf_mul_ui(fact, fact, i+1);

	// Do the first k-1 bits.  The last bit is always one.
	int mstart = mderiv-1;
	if (mstart < 0) mstart = 0;
	for (int i=mstart; i<order; i++)
	{
		if (digs[i])
		{
			cpx_times_mpf(term, zetan, fact);
			cpx_add(sum, sum, term);
		}
		cpx_mul(zetan, zetan, zeta);
		if (0 < mderiv)
		{
			mpf_mul_ui(fact, fact, i+2);
			mpf_div_ui(fact, fact, i-mderiv+2);
		}
	}
	// The final bit is always one.
	cpx_times_mpf(term, zetan, fact);
	cpx_add(sum, sum, term);

	// Subtract one.
	if (0 == mderiv)
		cpx_sub_ui(sum, sum, 1, 0);

	cpx_clear(zetan);
	cpx_clear(term);
	mpf_clear(fact);
}

// -------------------------------------------------------

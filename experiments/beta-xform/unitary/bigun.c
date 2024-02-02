/*
 * bigun.c
 *
 * Compute "unitary" zeros using bignum.
 * Feb 2024
 */

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

/**
 * Compute zeta polynomial. Roots of this poly will be the desired
 * zeros that we hope are accumulating.  This is the E(beta;z) thing.
 *
 * sum: returned value
 * zeta: location at which to evaluate.
 * digs: arry of digits to use
 * order: max order of the poly. There must be at least this many
 *     binary digits available.
 */
void ebz(cpx_t sum, cpx_t zeta, char* digs, int order)
{
	cpx_t zetan;
	cpx_init(zetan);
	cpx_set_ui(zetan, 1, 0);

	cpx_set_ui(sum, 0, 0);

	// Do the first k-1 of them.  The last one is always one.
	for (int i=0; i<order; i++)
	{
		if (digs[i]) cpx_add(sum, sum, zetan);
		cpx_mul(zetan, zetan, zeta);
	}
	// The final bit is always one.
	cpx_add(sum, sum, zetan);

	// times zeta
	cpx_mul(sum, sum, zeta);

	// subtract one.
	cpx_sub_ui(sum, sum, 1, 0);

	cpx_clear(zetan);
}

int main(int argc, char* argv[])
{
	int bprec = 500;
	mpf_set_default_prec(bprec);

	// Set beta to exactly 1.6
	mpf_t beta;
	mpf_init(beta);
	mpf_set_ui(beta, 16);
	mpf_div_ui(beta, beta, 10);

#define NBITS 400
	char bitseq[NBITS];
	gen_bitseq(beta, bitseq, NBITS);
	for (int i=0; i<70; i++)
		printf("%d", bitseq[i]);
	printf("\n");

	cpx_t zeta;
	cpx_init(zeta);
	cpx_set_ui(zeta, 1, 1);
	cpx_div_mpf(zeta, zeta, beta);

	cpx_t poly;
	cpx_init(poly);
	ebz(poly, zeta, bitseq, 100);
}
